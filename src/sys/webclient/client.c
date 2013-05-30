
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

#if defined(PETSC_USE_CERTIFICATE)
static int password_cb(char *buf,int num, int rwflag,void *userdata)
{
  if (num < strlen(PASSWORD)+1) return(0);
  strcpy(buf,PASSWORD);
  return(strlen(PASSWORD));
}
#endif

static void sigpipe_handle(int x)
{
}

#undef __FUNCT__
#define __FUNCT__ "PetscSSLInitializeContext"
/*
    PetscSSLInitializeContext - Set up an SSL context suitable for initiating HTTPS requests.

    If built with PETSC_USE_CERTIFICATE requires the user have created a self-signed certificate with

$    ./CA.pl  -newcert  (using the passphrase of password)
$    cat newkey.pem newcert.pem > sslclient.pem

    and put the resulting file in either the current directory (with the application) or in the home directory. This seems kind of
    silly but it was all I could figure out.

*/
PetscErrorCode PetscSSLInitializeContext(SSL_CTX **octx)
{
    SSL_METHOD     *meth;
    SSL_CTX        *ctx;
#if defined(PETSC_USE_CERTIFICATE)
    char           keyfile[PETSC_MAX_PATH_LEN];
    PetscBool      exists;
    PetscErrorCode ierr;
#endif

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

#if defined(PETSC_USE_CERTIFICATE)
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
#endif

    *octx = ctx;
    PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscSSLDestroyContext"
PetscErrorCode PetscSSLDestroyContext(SSL_CTX *ctx)
{
  PetscFunctionBegin;
  SSL_CTX_free(ctx);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "PetscHTTPSRequest"
/*
     PetscHTTPSRequest - Send a request to an HTTPS server

   Input Parameters:
+   type - either "POST" or "GET"
.   url - complete URL of request including https://
.   header - additional header information, may be NULL
.   ctype - data type of body, for example application/json
.   body - data to send to server
.   ssl - obtained with PetscHTTPSConnect()
-   buffsize - size of buffer

   Output Parameter:
.   buff - everything returned from server
 */
static PetscErrorCode PetscHTTPSRequest(const char type[],const char url[],const char header[],const char ctype[],const char body[],SSL *ssl,char buff[],size_t buffsize)
{
  char           *request=0;
  char           contentlength[40],contenttype[80];
  int            r;
  size_t         request_len,len,headlen,bodylen,contentlen,urllen,typelen,contenttypelen = 0;
  PetscErrorCode ierr;
  PetscBool      flg;

  PetscFunctionBegin;
  ierr = PetscStrbeginswith(url,"https://",&flg);CHKERRQ(ierr);
  if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"URL must begin with https://");
  if (header) {
    ierr = PetscStrendswith(header,"\r\n",&flg);CHKERRQ(ierr);
    if (!flg) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONG,"header must end with \\r\\n");
  }

  ierr = PetscStrlen(type,&typelen);CHKERRQ(ierr);
  ierr = PetscStrlen(url,&urllen);CHKERRQ(ierr);
  if (ctype) {
    ierr = PetscSNPrintf(contenttype,80,"Content-Type: %s\r\n",ctype);CHKERRQ(ierr);
    ierr = PetscStrlen(contenttype,&contenttypelen);CHKERRQ(ierr);
  }
  ierr = PetscStrlen(header,&headlen);CHKERRQ(ierr);
  ierr = PetscStrlen(body,&bodylen);CHKERRQ(ierr);
  ierr = PetscSNPrintf(contentlength,40,"Content-Length: %d\r\n\r\n",(int)bodylen);CHKERRQ(ierr);
  ierr = PetscStrlen(contentlength,&contentlen);CHKERRQ(ierr);

  /* Now construct our HTTP request */
  request_len = typelen + 1 + urllen + 35 + headlen + contenttypelen + contentlen + bodylen + 1;
  ierr = PetscMalloc(request_len*sizeof(char),&request);CHKERRQ(ierr);
  ierr = PetscStrcpy(request,type);CHKERRQ(ierr);
  ierr = PetscStrcat(request," ");CHKERRQ(ierr);
  ierr = PetscStrcat(request,url);CHKERRQ(ierr);
  ierr = PetscStrcat(request," HTTP/1.0\r\nUser-Agent:PETScClient\r\n");CHKERRQ(ierr);
  ierr = PetscStrcat(request,header);CHKERRQ(ierr);
  if (ctype) {
    ierr = PetscStrcat(request,contenttype);CHKERRQ(ierr);
  }
  ierr = PetscStrcat(request,contentlength);CHKERRQ(ierr);
  ierr = PetscStrcat(request,body);CHKERRQ(ierr);
  ierr = PetscStrlen(request,&request_len);CHKERRQ(ierr);
  ierr = PetscInfo1(NULL,"HTTPS request follows: \n%s\n",request);CHKERRQ(ierr);

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

#undef __FUNCT__
#define __FUNCT__ "PetscHTTPSConnect"
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

#undef __FUNCT__
#define __FUNCT__ "PetscURLShorten"
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
  ierr = PetscHTTPSRequest("POST","https://www.googleapis.com/urlshortener/v1/url",NULL,"application/json",body,ssl,buff,1024);CHKERRQ(ierr);
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


