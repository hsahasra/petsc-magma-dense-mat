
#include <petscsys.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <netdb.h>
#include <fcntl.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

#include <openssl/ssl.h>
#include <openssl/err.h>

static BIO *bio_err = NULL;

#define PASSWORD "password"

static int password_cb(char *buf,int num, int rwflag,void *userdata)
{
  if (num < strlen(PASSWORD)+1) return(0);
  strcpy(buf,PASSWORD);
  return(strlen(PASSWORD));
}

static void sigpipe_handle(int x)
{
}

/*
    PetscSSLInitializeContext - Set up an SSL context suitable for initiating HTTPS requests.

    Requires the user have created a self-signed certificate with

$    ./CA.pl  -newcert  (using the passphrase of password)
$    cat newkey.pem newcert.pem > sslclient.pem

    and put the resulting file in either the current directory (with the application) or in the home directory. This seems kind of
    silly but it was all I could figure out.

*/
PetscErrorCode PetscSSLInitializeContext(SSL_CTX **octx)
{
    SSL_METHOD     *meth;
    SSL_CTX        *ctx;
    char           keyfile[PETSC_MAX_PATH_LEN];
    PetscBool      exists;
    PetscErrorCode ierr;

    PetscFunctionBegin;
    if (!bio_err){
      SSL_library_init();
      SSL_load_error_strings();
      bio_err = BIO_new_fp(stderr,BIO_NOCLOSE);
    }

    /* Set up a SIGPIPE handler */
    signal(SIGPIPE,sigpipe_handle);

    meth = SSLv23_method();
    ctx  = SSL_CTX_new(meth);

    /* Locate keyfile */
    ierr = PetscStrcpy(keyfile,"sslclient.pem");CHKERRQ(ierr);
    ierr = PetscTestFile(keyfile,'r',&exists);CHKERRQ(ierr);
    if (!exists) {
      ierr = PetscGetHomeDirectory(keyfile,PETSC_MAX_PATH_LEN);CHKERRQ(ierr);
      ierr = PetscStrcat(keyfile,"/");CHKERRQ(ierr);
      ierr = PetscStrcat(keyfile,"sslclient.pem");CHKERRQ(ierr);
      ierr = PetscTestFile(keyfile,'r',&exists);CHKERRQ(ierr);
      if (!exists) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Unable to locate sslclient.pem file in current directory or home directory");
    }

    /* Load our keys and certificates*/
    if (!(SSL_CTX_use_certificate_chain_file(ctx,keyfile))) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot read certificate file");

    SSL_CTX_set_default_passwd_cb(ctx,password_cb);
    if (!(SSL_CTX_use_PrivateKey_file(ctx,keyfile,SSL_FILETYPE_PEM))) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_FILE_OPEN,"Cannot read key file");

    *octx = ctx;
    PetscFunctionReturn(0);
}

PetscErrorCode PetscSSLDestroyContext(SSL_CTX *ctx)
{
  PetscFunctionBegin;
  SSL_CTX_free(ctx);
  PetscFunctionReturn(0);
}

static PetscErrorCode PetscHTTPSRequest(const char header[],const char body[],SSL *ssl,char buff[],size_t buffsize)
{
  char           *request=0;
  char           contentlength[40];
  int            r;
  size_t         request_len,len,headlen,bodylen,contentlen;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = PetscStrlen(header,&headlen);CHKERRQ(ierr);
  ierr = PetscStrlen(body,&bodylen);CHKERRQ(ierr);
  ierr = PetscSNPrintf(contentlength,40,"Content-Length: %d\r\n\r\n",(int)bodylen);CHKERRQ(ierr);
  ierr = PetscStrlen(contentlength,&contentlen);CHKERRQ(ierr);

  /* Now construct our HTTP request */
  request_len = headlen + bodylen + contentlen + 1;
  ierr = PetscMalloc(request_len*sizeof(char),&request);CHKERRQ(ierr);
  ierr = PetscStrcpy(request,header);CHKERRQ(ierr);
  ierr = PetscStrcat(request,contentlength);CHKERRQ(ierr);
  ierr = PetscStrcat(request,body);CHKERRQ(ierr);
  ierr = PetscStrlen(request,&request_len);CHKERRQ(ierr);
  ierr = PetscInfo1(NULL,"HTTPS request: %s\n",request);CHKERRQ(ierr);

  r = SSL_write(ssl,request,request_len);
  switch (SSL_get_error(ssl,r)){
    case SSL_ERROR_NONE:
      if (request_len != r) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Incomplete write to SSL socket");
      break;
      default:
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"SSL socket write problem");
  }

  /* Now read the server's response, assuming  that it's terminated by a close */
  while (1){
    r = SSL_read(ssl,buff,(int)buffsize);
    switch (SSL_get_error(ssl,r)){
      case SSL_ERROR_NONE:
        len = r;
        break;
      case SSL_ERROR_ZERO_RETURN:
        goto shutdown;
      case SSL_ERROR_SYSCALL:
        goto done;
      default:
        SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"SSL read problem");
    }
    buff[len] = 0; /* null terminate string */
  }

  shutdown:
    SSL_shutdown(ssl);  /* ignore shutdown error message */

  done:
    SSL_free(ssl);
    ierr = PetscFree(request);CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

PetscErrorCode PetscHTTPSConnect(const char host[],int port,SSL_CTX *ctx,int *sock,SSL **ssl)
{
  BIO            *sbio;
  PetscErrorCode ierr;

  PetscFunctionBegin;
  /* Connect the TCP socket*/
  ierr = PetscOpenSocket(host,port,sock);CHKERRQ(ierr);

  /* Connect the SSL socket */
  *ssl = SSL_new(ctx);
  sbio = BIO_new_socket(*sock,BIO_NOCLOSE);
  SSL_set_bio(*ssl,sbio,sbio);
  if (SSL_connect(*ssl) <= 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"SSL connect error");
  PetscFunctionReturn(0);
}

/*@C
     PetscURLShorten - Uses Google's service to get a short url for a long url

    Input Parameters:
+    url - long URL you want shortened
-    lenshorturl - length of buffer to contain short URL

    Output Parameter:
.    shorturl - the shortened URL

@*/
PetscErrorCode PetscURLShorten(const char url[],char shorturl[],size_t lenshorturl)
{
  SSL_CTX        *ctx;
  SSL            *ssl;
  int            sock;
  PetscErrorCode ierr;
  char           buff[1024],body[512],*sub1,*sub2;

  PetscFunctionBegin;
  ierr = PetscSSLInitializeContext(&ctx);CHKERRQ(ierr);
  ierr = PetscHTTPSConnect("www.googleapis.com",443,ctx,&sock,&ssl);CHKERRQ(ierr);
  ierr = PetscSNPrintf(body,512,"{\"longUrl\": \"%s\"}",url);CHKERRQ(ierr);
  ierr = PetscHTTPSRequest("POST https://www.googleapis.com/urlshortener/v1/url HTTP/1.0\r\nUser-Agent:PETScClient\r\nContent-type: application/json\r\n",body,ssl,buff,1024);CHKERRQ(ierr);
  PetscSSLDestroyContext(ctx);
  close(sock);
  ierr = PetscInfo1(NULL,"Response from Google shortener %s\n",buff);
  ierr = PetscStrstr(buff,"\"id\": \"",&sub1);CHKERRQ(ierr);
  if (sub1) {
    sub1 += 7;
    ierr = PetscStrstr(sub1,"\"",&sub2);CHKERRQ(ierr);
    if (!sub2) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Google did not shorten URL");
    sub2[0] = 0;
    ierr = PetscStrncpy(shorturl,sub1,lenshorturl);CHKERRQ(ierr);
  } else SETERRQ(PETSC_COMM_SELF,PETSC_ERR_LIB,"Google did not shorten URL");
  PetscFunctionReturn(0);
}


