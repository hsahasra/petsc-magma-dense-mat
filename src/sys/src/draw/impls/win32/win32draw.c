/* $Id: win32draw.c,v 1.5 2000/08/03 20:03:15 balay Exp balay $ */
#include "petsc.h"
#include "src/sys/src/draw/drawimpl.h"
#include "win32draw.h"

#define IDC_FOUR       109
#define IDI_FOUR       107
#define IDM_EXIT       105
#define IDR_POPUP      103
#define MAX_LOADSTRING 100

#define SelectPen(hdc, hpen)  ((HPEN)SelectObject((hdc),  (HGDIOBJ)(HPEN)(hpen)))
#define SelectFont(hdc,hfont) ((HFONT)SelectObject((hdc), (HGDIOBJ)(HFONT)(hfont)))
#define GetStockBrush(i)      ((HBRUSH)GetStockObject(i))

#define XTRANS(draw,win,x) \
   (int)(((win)->w)*((draw)->port_xl + (((x - (draw)->coor_xl)*\
                                   ((draw)->port_xr - (draw)->port_xl))/\
                                   ((draw)->coor_xr - (draw)->coor_xl))))
#define YTRANS(draw,win,y) \
   (int)(((win)->h)*(1.0-(draw)->port_yl - (((y - (draw)->coor_yl)*\
                                   ((draw)->port_yr - (draw)->port_yl))/\
                                   ((draw)->coor_yr - (draw)->coor_yl))))

HINSTANCE     hInst;
HANDLE        g_hWindowListMutex = NULL;
WindowNode    WindowListHead     = NULL;

/* Hard coded color hue until hue.c works with this */
unsigned char RedMap[]   = {255,0,255,0,0,0,255,127,34,255,238,165,255,255,190,255,255,238,0,255,105,154,135,0,0,244,152,176,220,216,50,255};
unsigned char GreenMap[] = {255,0,0,255,255,0,0,255,139,165,130,42,182,127,190,255,215,162,197,246,105,205,206,100,0,164,245,224,17,191,205,240};
unsigned char BlueMap[]  = {255,0,0,0,255,255,225,212,34,0,238,42,193,80,190,0,0,173,205,143,105,50,235,0,128,96,255,230,120,216,50,245};

/* Foward declarations of functions included in this code module: */
LRESULT CALLBACK  WndProc(HWND, UINT, WPARAM, LPARAM);
extern int  TranslateColor_Win32(Draw,int);
extern int  AverageColorRectangle_Win32(Draw,int,int,int,int);
extern int  AverageColorTriangle_Win32(Draw,int,int,int);
extern void MessageLoopThread(Draw_Win32 *);
extern int  deletemouselist_Win32(WindowNode);
extern void OnPaint_Win32(HWND);
extern void OnDestroy_Win32(HWND);
extern void OnSize_Win32(HWND,UINT,int,int);
extern int  MouseRecord_Win32(HWND,DrawButton);
extern int  DrawGetPopup_Win32(Draw,Draw *);

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawSetDoubleBuffer_Win32"></a> */"DrawSetDoubleBuffer_Win32" 
static int DrawSetDoubleBuffer_Win32(Draw draw)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HDC        hdc      = GetDC(windraw->hWnd);
  
  PetscFunctionBegin;
  windraw->node->DoubleBuffer = CreateCompatibleDC(hdc);
  windraw->node->DoubleBufferBit = CreateCompatibleBitmap(hdc,windraw->w,windraw->h);
  windraw->node->dbstore = SelectObject(windraw->node->DoubleBuffer,windraw->node->DoubleBufferBit);
  /* Fill background of second buffer */
  ExtFloodFill(windraw->node->DoubleBuffer,0,0,COLOR_WINDOW,FLOODFILLBORDER);
  /* Copy current buffer into seconf buffer and set window data as double buffered */
  BitBlt(windraw->node->DoubleBuffer,
         0,0,
         windraw->w,windraw->h,
         windraw->node->Buffer,
         0,0,
         SRCCOPY);

  windraw->node->DoubleBuffered = PETSC_TRUE;
  ReleaseDC(windraw->hWnd,hdc);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawFlush_Win32"></a> */"DrawFlush_Win32" 
static int DrawFlush_Win32(Draw draw)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HDC        hdc = GetDC(windraw->hWnd);
  
  PetscFunctionBegin;
  /* flush double buffer into primary buffer */
  BitBlt(windraw->node->Buffer,
         0,0,
         windraw->w,windraw->h,
         windraw->node->DoubleBuffer,
         0,0,
         SRCCOPY);
  /* flush double buffer into window */
  BitBlt(hdc,
         0,0,
         windraw->w,windraw->h,
         windraw->node->DoubleBuffer,
         0,0,
         SRCCOPY);
  ReleaseDC(windraw->hWnd,hdc);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="deletemouselist_Win32"></a> */"deletemouselist_Win32" 
static int deletemouselist_Win32(WindowNode deletelist)
{ 
  /* Called upon window close. Frees memory of linked list of stored mouse commands */
  MouseNode node;
  
  while(deletelist->MouseListHead != NULL) {       
    node = deletelist->MouseListHead;
    if(deletelist->MouseListHead->mnext != NULL) {
      deletelist->MouseListHead = deletelist->MouseListHead->mnext;
    }
    PetscFree(node);
  }
  deletelist->MouseListHead = deletelist->MouseListTail = NULL;
  if (deletelist->wprev != NULL) {
    deletelist->wprev->wnext = deletelist->wnext;
  }
  if (deletelist->wnext != NULL) {
    deletelist->wnext->wprev = deletelist->wprev;
  }
  PetscFree(deletelist);
  return 0;
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawGetMouseButton_Win32"></a> */"DrawGetMouseButton_Win32" 
static int DrawGetMouseButton_Win32(Draw draw, DrawButton *button,double *x_user,double *y_user, \
                                    double *x_phys,double *y_phys)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  WindowNode  current;
  MouseNode   node=0;
  
  PetscFunctionBegin;
  /* Make sure no other code is using the linked list at this moment */
  WaitForSingleObject(g_hWindowListMutex, INFINITE);
  /* Look for the node that matches the window you are using */
  current = WindowListHead;
  while (current != NULL) {
    if(current->hWnd == windraw->hWnd) {       
      current->IsGetMouseOn = TRUE;
      break;
    } else {
      current = current->wnext;
    }
  }
  /* If no actions have occured, wait for one */
  node = current->MouseListHead;
  if (node == NULL) {
    ReleaseMutex(g_hWindowListMutex);
    WaitForSingleObject(current->event, INFINITE);
    WaitForSingleObject(g_hWindowListMutex, INFINITE);
  }
  /* once we have the information, assign the pointers to it */
  *button = current->MouseListHead->Button;
  *x_user = current->MouseListHead->user.x;
  *y_user = current->MouseListHead->user.y;
  /* optional arguments */
  if (x_phys) *x_phys = current->MouseListHead->phys.x;
  if (y_phys) *y_phys = current->MouseListHead->phys.y;
  /* remove set of information from sub linked-list, delete the node */
  current->MouseListHead = current->MouseListHead->mnext;
  if (current->MouseListHead == NULL) {
    ResetEvent(current->event);
    current->MouseListTail = NULL;
  }
  if (node) PetscFree(node);

  /* Release mutex so that  other code can use
     the linked list now that we are done with it */
  ReleaseMutex(g_hWindowListMutex);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawPause_Win32"></a> */"DrawPause_Win32" 
static int DrawPause_Win32(Draw draw)
{
  PetscFunctionBegin;
  PetscSleep(draw->pause);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="TranslateColor_Win32"></a> */"TranslateColor_Win32" 
static int TranslateColor_Win32(Draw draw,int color)
{
  /* Maps single color value into the RGB colors in our tables */
  Draw_Win32 *windraw   = (Draw_Win32*)draw->data;
  windraw->currentcolor = RGB(RedMap[color],GreenMap[color],BlueMap[color]);
  return 0;
}
#undef __FUNC__  
#define __FUNC__ /* <a name="AverageColorRectangle_Win32"></a> */"AverageColorRectangle_Win32"
static int AverageColorRectangle_Win32(Draw draw,int c1,int c2, int c3, int c4)
{
  /* Averages colors given at points of rectangle and sets color from color table
    will be changed once the color gradient problem is worked out */
  Draw_Win32 *windraw   = (Draw_Win32*)draw->data;
  windraw->currentcolor = RGB(((RedMap[c1]+RedMap[c2]+RedMap[c3]+RedMap[c4])/4),
                              ((GreenMap[c1]+GreenMap[c2]+GreenMap[c3]+GreenMap[c4])/4),
                              ((BlueMap[c1]+BlueMap[c2]+BlueMap[c3]+BlueMap[c4])/4));
  return 0;
}
#undef __FUNC__  
#define __FUNC__ /* <a name="AverageColorTriangle_Win32"></a> */"AverageColorTriangle_Win32"
static int AverageColorTriangle_Win32(Draw draw,int c1,int c2,int c3)
{
  /* Averages colors given at points of rectangle and sets color from color table
    will be changed once the color gradient problem is worked out */
  Draw_Win32 *windraw   = (Draw_Win32*)draw->data;
  windraw->currentcolor = RGB((RedMap[c1]+RedMap[c2]+RedMap[c3])/3,
                              (GreenMap[c1]+GreenMap[c2]+GreenMap[c3])/3,
                              (BlueMap[c1]+BlueMap[c2]+BlueMap[c3])/3); 
  return 0;
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawRectangle_Win32"></a> */"DrawRectangle_Win32"
static int DrawRectangle_Win32(Draw draw,double xl,double yl,double xr,double yr,int c1,int c2,int c3,int c4)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HBRUSH      hbrush;
  RECT        rect;
  int         x1,y1,x2,y2;
  HDC         hdc;
  
  PetscFunctionBegin;
  x1 = XTRANS(draw,windraw,xl);
  x2 = XTRANS(draw,windraw,xr);
  y1 = YTRANS(draw,windraw,yl);
  y2 = YTRANS(draw,windraw,yr);
  SetRect(&rect,x1,y2,x2,y1);        
  if (c1==c2 && c2==c3 && c3==c4) {         
    TranslateColor_Win32(draw,c1);
  } else {                                   
    AverageColorRectangle_Win32(draw,c1,c2,c3,c4);
  }
  hbrush = CreateSolidBrush(windraw->currentcolor);
  
  if(windraw->node->DoubleBuffered) {
    hdc = windraw->node->DoubleBuffer;
  } else {
    hdc = windraw->node->Buffer;
  }
  FillRect(hdc,&rect,hbrush);
  /* Forces a WM_PAINT message and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawLine_Win32"></a> */"DrawLine_Win32"
static int DrawLine_Win32(Draw draw,double xl,double yl,double xr,double yr,int color)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HPEN       hpen;
  int        x1,y1,x2,y2;
  HDC        hdc;
  
  PetscFunctionBegin;
  TranslateColor_Win32(draw,color);
  x1   = XTRANS(draw,windraw,xl);x2  = XTRANS(draw,windraw,xr); 
  y1   = YTRANS(draw,windraw,yl);y2  = YTRANS(draw,windraw,yr); 
  hpen = CreatePen (PS_SOLID, windraw->linewidth, windraw->currentcolor);
  if(windraw->node->DoubleBuffered) {
    hdc = windraw->node->DoubleBuffer;
  } else {
    hdc = windraw->node->Buffer;
  }
  SelectPen(hdc,hpen);
  MoveToEx(hdc,x1,y1,NULL);
  LineTo(hdc,x2,y2);
  /* Forces a WM_PAINT message and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawLineSetWidth_Win32"></a> */"DrawLineSetWidth_Win32"
static int DrawLineSetWidth_Win32(Draw draw,double width)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  int         averagesize,finalwidth;
  RECT        rect;
  
  PetscFunctionBegin;
  GetClientRect(windraw->hWnd,&rect);
  averagesize = ((rect.right - rect.left)+(rect.bottom - rect.top))/2;
  finalwidth  = (int)floor(averagesize*width);
  if (finalwidth < 1) {
    finalwidth = 1; /* minimum size DrawLine can except */
  }
  windraw->linewidth = finalwidth;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawLineGetWidth_Win32"></a> */"DrawLineGetWidth_Win32"
static int DrawLineGetWidth_Win32(Draw draw,PetscReal *width)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  *width = (PetscReal)windraw->linewidth;
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawPoint_Win32"></a> */"DrawPoint_Win32"
static int DrawPoint_Win32(Draw draw,double x,double y,int color)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HBRUSH     hbrush;
  HRGN       hrgn;
  int        radius;
  int        x1,y1,left,right,top,bottom;
  HDC        hdc;
  
  PetscFunctionBegin;
  TranslateColor_Win32(draw,color);
  x1     = XTRANS(draw,windraw,x);   
  y1     = YTRANS(draw,windraw,y);
  hbrush = CreateSolidBrush(windraw->currentcolor);
  if(windraw->node->DoubleBuffered) {
    hdc = windraw->node->DoubleBuffer;
  } else {
    hdc = windraw->node->Buffer;
  }
  /* desired size is one logical pixel so just turn it on */
  if (windraw->pointdiameter == 1) {
    SetPixelV(hdc,x1,y1,windraw->currentcolor);
  } else {
    /* draw point around position determined */
    radius = (int)floor(windraw->pointdiameter/2);
    hrgn   = CreateEllipticRgn(x1-radius,y1-radius,x1+radius,y1+radius);
    FillRgn(hdc,hrgn,hbrush);
  }
  /* Forces a WM_PAINT and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawPointSetSize_Win32"></a> */"DrawPointSetSize_Win32"
static int DrawPointSetSize_Win32(Draw draw,double width)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  int         averagesize,diameter;
  RECT        rect;
  
  PetscFunctionBegin;
  GetClientRect(windraw->hWnd,&rect);
  averagesize = ((rect.right - rect.left)+(rect.bottom - rect.top))/2;
  diameter    = (int)floor(averagesize*width);
  if (diameter < 1) diameter = 1;
  windraw->pointdiameter     = diameter;
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawString_Win32"></a> */"DrawString_Win32"
static int DrawString_Win32(Draw draw,double x,double y,int color,char *text)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  RECT        r;
  HFONT       hfont;                                                                    
  LOGFONT     logfont; 
  int         x1,y1;
  HDC         hdc;
  
  PetscFunctionBegin;
  x1              = XTRANS(draw,windraw,x);
  y1              = YTRANS(draw,windraw,y);
  r.bottom        = y1;
  r.left          = x1;
  r.right         = x1 + 1; 
  r.top           = y1 + 1;
  logfont.lfHeight         = windraw->stringheight;
  logfont.lfWidth          = windraw->stringwidth;
  logfont.lfEscapement     = 0;
  logfont.lfOrientation    = 0;
  logfont.lfCharSet        = 0;
  logfont.lfClipPrecision  = 0;
  logfont.lfItalic         = 0;
  logfont.lfOutPrecision   = 0;
  logfont.lfPitchAndFamily = DEFAULT_PITCH;
  logfont.lfQuality        = DEFAULT_QUALITY;
  logfont.lfStrikeOut      = 0;
  logfont.lfUnderline      = 0;
  logfont.lfWeight         = FW_NORMAL;
  hfont = CreateFontIndirect(&logfont); 
  TranslateColor_Win32(draw,color);
  if(windraw->node->DoubleBuffered) {
    hdc = windraw->node->DoubleBuffer;
  } else {
    hdc = windraw->node->Buffer;
  }
  SelectFont(hdc,hfont);
  SetTextColor(hdc,windraw->currentcolor);
  DrawText(hdc,text,lstrlen(text),&r,DT_NOCLIP);
  DeleteObject(hfont);
  /* Forces a WM_PAINT message and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawStringVertical_Win32"></a> */"DrawStringVertical_Win32"
static int DrawStringVertical_Win32(Draw draw,double x,double y,int color,char *text)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  RECT        r;
  HFONT       hfont;                                                                                    
  LOGFONT     logfont;
  int         x1,y1;
  HDC         hdc;
  
  PetscFunctionBegin;
  x1           = XTRANS(draw,windraw,x);
  y1           = XTRANS(draw,windraw,y);
  r.bottom     = x1;
  r.left       = y1 + 30;
  r.right      = x1 + 1;
  r.top        = y1 - 30;
  logfont.lfEscapement     = 2700; /* Causes verticle text drawing */
  logfont.lfHeight         = windraw->stringheight;
  logfont.lfWidth          = windraw->stringwidth;
  logfont.lfOrientation    = 0;
  logfont.lfCharSet        = DEFAULT_CHARSET;
  logfont.lfClipPrecision  = 0;
  logfont.lfItalic         = 0;
  logfont.lfOutPrecision   = 0;
  logfont.lfPitchAndFamily = DEFAULT_PITCH;
  logfont.lfQuality        = DEFAULT_QUALITY;
  logfont.lfStrikeOut      = 0;
  logfont.lfUnderline      = 0;
  logfont.lfWeight         = FW_NORMAL;
  hfont = CreateFontIndirect(&logfont);
  TranslateColor_Win32(draw,color);
  if(windraw->node->DoubleBuffered) {
    hdc = windraw->node->DoubleBuffer;
  } else {
    hdc = windraw->node->Buffer;
  }
  SelectFont(hdc,hfont);
  SetTextColor(hdc,windraw->currentcolor);
  DrawText(hdc,text,lstrlen(text),&r,DT_NOCLIP | DT_SINGLELINE );
  DeleteObject(hfont);
  /* Forces a WM_PAINT message and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawStringSetSize_Win32"></a> */"DrawStringSetSize_Win32"
static int DrawStringSetSize_Win32(Draw draw,double width,double height)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  int         w,h;
  
  PetscFunctionBegin;
  w = (int)((windraw->w)*width *(draw->port_xr - draw->port_xl)/(draw->coor_xr - draw->coor_xl));
  h = (int)((windraw->h)*height*(draw->port_yr - draw->port_yl)/(draw->coor_yr - draw->coor_yl));
  if (h < 1) h = 1;
  if (w < 1) w = 1;
  windraw->stringheight = h;
  windraw->stringwidth  = w;
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* a name="DrawStringGetSize_Win32"></a> */"DrawStringGetSize_Win32"
static int DrawStringGetSize_Win32(Draw draw,double *width,double *height)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  *height = (double)windraw->stringheight;
  *width  = (double)windraw->stringwidth;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawResizeWindow_Win32"></a> */"DrawResizeWindow_Win32"
static int DrawResizeWindow_Win32(Draw draw,int w,int h)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  RECT        r;
  
  PetscFunctionBegin;
  GetWindowRect(windraw->hWnd,&r);
  MoveWindow(windraw->hWnd,r.left,r.top,(int)w,(int)h,TRUE);
  /* set all variable dealing with window dimensions */
  windraw->node->bitheight = windraw->h = draw->h = h;
  windraw->node->bitwidth  = windraw->w = draw->w = w;
  /* set up graphic buffers with the new size of window */
  SetBitmapDimensionEx(windraw->node->BufferBit,w,h,NULL);
  if(windraw->node->DoubleBuffered) {
    SetBitmapDimensionEx(windraw->node->DoubleBufferBit,w,h,NULL);
  }
  windraw->haveresized = PETSC_TRUE;
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawCheckResizeWindow_Win32"></a> */"DrawCheckResizeWindow_Win32"
static int DrawCheckResizedWindow_Win32(Draw draw)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  if (windraw->haveresized == 1) {
    PetscFunctionReturn(1);
  } else {
    PetscFunctionReturn(0);
  }
  
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawSetTitle_Win32"></a> */"DrawSetTitle_Win32"
static int DrawSetTitle_Win32(Draw draw, char *title)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  draw->title = title;
  SetWindowText(windraw->hWnd,title);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawClear_Win32"></a> */"DrawClear_Win32"
static int DrawClear_Win32(Draw draw)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  /* clear primary buffer */
  ExtFloodFill(windraw->node->Buffer,0,0,COLOR_WINDOW,FLOODFILLBORDER);
  /* if exists clear secondary buffer */
  if(windraw->node->DoubleBuffered) {
    ExtFloodFill(windraw->node->DoubleBuffer,0,0,COLOR_WINDOW,FLOODFILLBORDER);
  }
  /* force WM_PAINT message so cleared buffer will show */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawTriangle_Win32"></a> */"DrawTriangle_Win32"
static int DrawTriangle_Win32(Draw draw,double x1,double y1,double x2,double y2,double x3,double y3,
                              int c1,int c2,int c3)
{       
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  HBRUSH      hbrush;
  HPEN        hpen;
  int         p1x,p1y,p2x,p2y,p3x,p3y;
  HDC         bit;
  
  PetscFunctionBegin;
  AverageColorTriangle_Win32(draw,c1,c2,c3); 
  hbrush = CreateSolidBrush(windraw->currentcolor);
  hpen   = CreatePen(PS_SOLID,0,windraw->currentcolor);
  p1x = XTRANS(draw,windraw,x1);
  p2x = XTRANS(draw,windraw,x2);
  p3x = XTRANS(draw,windraw,x3);
  p1y = YTRANS(draw,windraw,y1);
  p2y = YTRANS(draw,windraw,y2);
  p3y = YTRANS(draw,windraw,y3);
  
  if(windraw->node->DoubleBuffered) {
    bit = windraw->node->DoubleBuffer;
  } else {
    bit = windraw->node->Buffer;
  }
  BeginPath(bit);
  MoveToEx(bit,p1x,p1y,NULL);
  LineTo(bit,p2x,p2y);
  LineTo(bit,p3x,p3y);
  LineTo(bit,p1x,p1y);
  EndPath(bit);
  SelectPen(bit,hpen);
  SelectBrush(bit,hbrush);
  StrokeAndFillPath(bit);
  /* Forces a WM_PAINT message and erases background */
  InvalidateRect(windraw->hWnd,NULL,TRUE);
  UpdateWindow(windraw->hWnd);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="PopMessageLoopThread_Win32"></a> */"PopMessageLoopThread_Win32"
void PopMessageLoopThread_Win32(Draw popdraw)
{
  Draw_Win32 *pop = (Draw_Win32*)popdraw->data;
  MSG         msg;
  HWND        hWnd = NULL;
  char        PopClassName [MAX_LOADSTRING + 1]; 
  RECT        r;
  int         width,height;
  WNDCLASSEX  myclass;
  LPVOID      lpMsgBuf;
  
  /* initialize window class parameters */
  myclass.cbSize        = sizeof(WNDCLASSEX);
  myclass.style         = CS_OWNDC;
  myclass.lpfnWndProc   = (WNDPROC)WndProc;
  myclass.cbClsExtra    = 0;
  myclass.cbWndExtra    = 0;
  myclass.hInstance     = NULL;
  myclass.hIcon         = NULL;
  myclass.hCursor       = LoadCursor(NULL, IDC_ARROW);
  myclass.hbrBackground = (HBRUSH)(COLOR_WINDOW + 1);
  myclass.lpszMenuName  = NULL;
  myclass.lpszClassName = PopClassName;
  myclass.hIconSm       = NULL;
  
  RegisterClassEx(&myclass);
  
  SetRect(&r,0,0,450,450);
  
  width    = (r.right - r.left) / 3;
  height   = (r.bottom - r.top) / 3;
  
  hWnd = CreateWindowEx(0,
                        PopClassName,
                        NULL, 
                        WS_POPUPWINDOW | WS_CAPTION,
                        0,0, 
                        width,height,
                        NULL,
                        NULL,
                        hInst,
                        NULL);
  pop->x = 0;
  pop->y = 0;
  pop->w = width;
  pop->h = height;
  
  if(hWnd == NULL) {
    lpMsgBuf = "Window Not Succesfully Created";
    MessageBox( NULL, (LPCTSTR)lpMsgBuf, "Error", MB_OK | MB_ICONINFORMATION );
    LocalFree( lpMsgBuf );
    exit(0);
  }
  pop->hWnd = hWnd;
  /* display and update new popup window */
  ShowWindow(pop->hWnd, SW_SHOWNORMAL);
  UpdateWindow(pop->hWnd);
  SetEvent(pop->hReadyEvent);
  
  while (GetMessage(&msg, pop->hWnd, 0, 0)) {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
}

#undef __FUNC__  
#define __FUNC__ /* <a name="DrawDestroy_Win32"></a> */"DrawDestroy_Win32"
static int DrawDestroy_Win32(Draw draw)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  
  PetscFunctionBegin;
  SendMessage(windraw->hWnd,WM_DESTROY,0,0);
  PetscFree(windraw);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawSynchronizedFlush_Win32"></a> */"DrawSynchronizedFlush_Win32"
static int DrawSynchronizedFlush_Win32(Draw draw)
{
  /* Multi Processor is not implemeted yet */
  PetscFunctionBegin;
  DrawFlush_Win32(draw);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="DrawSynchronizedClear_Win32"></a> */"DrawSynchronizedClear_Win32"
static int DrawSynchronizedClear_Win32(Draw draw)
{
  /* Multi Processor is not implemeted yet */
  PetscFunctionBegin;
  DrawClear_Win32(draw);
  PetscFunctionReturn(0);
}
#undef __FUNC__  
#define __FUNC__ /* <a name="MessageLoopThread_Win32"></a> */"MessageLoopThread_Win32"
void MessageLoopThread_Win32(Draw draw)
{
  Draw_Win32 *windraw = (Draw_Win32*)draw->data;
  MSG        msg;
  HWND       hWnd = NULL;
  char       classname[MAX_LOADSTRING + 1];
  WNDCLASSEX wclass;
  LPVOID     lpMsgBuf;
  
  /* initialize window class parameters */
  wclass.cbSize         = sizeof(WNDCLASSEX);
  wclass.style          = CS_SAVEBITS | CS_HREDRAW | CS_VREDRAW;
  wclass.lpfnWndProc    = (WNDPROC)WndProc;
  wclass.cbClsExtra     = 0;
  wclass.cbWndExtra     = 0;
  wclass.hInstance      = NULL;
  wclass.hIcon          = LoadIcon(NULL,IDI_APPLICATION);
  wclass.hCursor        = LoadCursor(NULL,IDC_ARROW);
  wclass.hbrBackground  = GetStockBrush(WHITE_BRUSH);
  wclass.lpszMenuName   = NULL;
  wclass.lpszClassName  = classname;
  wclass.hIconSm        = NULL;
  
  RegisterClassEx(&wclass);
  
  
  hWnd = CreateWindowEx(0,
                        classname,
                        NULL,
                        WS_OVERLAPPEDWINDOW,
                        draw->x,
                        draw->y, 
                        draw->w,
                        draw->h, 
                        NULL,
                        NULL,
                        hInst,
                        NULL);
  
  if(hWnd == NULL) {
    lpMsgBuf = "Window Not Succesfully Created";
    MessageBox( NULL, (LPCTSTR)lpMsgBuf, "Error", MB_OK | MB_ICONINFORMATION );
    LocalFree( lpMsgBuf );
    exit(0);
  }
  windraw->hWnd = hWnd;
  /* display and update new window */
  ShowWindow(hWnd,SW_SHOWNORMAL);
  UpdateWindow(hWnd);
  SetEvent(windraw->hReadyEvent);
  
  while (GetMessage(&msg,hWnd, 0, 0)) {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }
}


static struct _DrawOps DvOps = { DrawSetDoubleBuffer_Win32,
                                 DrawFlush_Win32,
                                 DrawLine_Win32,
                                 DrawLineSetWidth_Win32,
                                 DrawLineGetWidth_Win32,
                                 DrawPoint_Win32,
                                 DrawPointSetSize_Win32,
                                 DrawString_Win32,
                                 DrawStringVertical_Win32,
                                 DrawStringSetSize_Win32,
                                 DrawStringGetSize_Win32,
                                 0,
                                 DrawClear_Win32,
                                 DrawSynchronizedFlush_Win32,
                                 DrawRectangle_Win32,
                                 DrawTriangle_Win32,
                                 DrawGetMouseButton_Win32,
                                 DrawPause_Win32,
                                 DrawSynchronizedClear_Win32,
                                 0,
                                 0,
                                 DrawGetPopup_Win32,
                                 DrawSetTitle_Win32,
                                 DrawCheckResizedWindow_Win32,
                                 DrawResizeWindow_Win32,
                                 DrawDestroy_Win32,
                                 0,
                                 0,
                                 0,
                                 0};


static int DrawGetPopup_Win32(Draw draw,Draw *popdraw)
{
  Draw_Win32 *pop;
  HANDLE      hThread = NULL;
  WindowNode  newnode;
  int         ierr;
  
  PetscFunctionBegin;
  pop = (Draw_Win32 *) PetscMalloc(sizeof (Draw_Win32));
  (*popdraw)->data = pop;
  
  /* the following is temporary fix for initializing a global datastructure */
  if(!g_hWindowListMutex) {
    g_hWindowListMutex = CreateMutex(NULL,FALSE,NULL);
  }
  ierr = PetscMemcpy((*popdraw)->ops,&DvOps,sizeof(DvOps));CHKERRQ(ierr);
  
  pop->hReadyEvent = CreateEvent(NULL, TRUE, FALSE, NULL);
  CreateThread(NULL, 0,(LPTHREAD_START_ROUTINE)PopMessageLoopThread_Win32,*popdraw,0,(unsigned long*)hThread);
  CloseHandle(hThread);
  WaitForSingleObject(pop->hReadyEvent, INFINITE);
  CloseHandle(pop->hReadyEvent);
  WaitForSingleObject(g_hWindowListMutex, INFINITE);
  
  draw->popup             = (*popdraw);
  newnode                 = (WindowNode)PetscMalloc(sizeof(struct _p_WindowNode));CHKPTRQ(newnode);
  newnode->MouseListHead  = NULL;
  newnode->MouseListTail  = NULL;
  newnode->wnext          = WindowListHead;
  newnode->wprev          = NULL;
  newnode->hWnd           = pop->hWnd;
  if(WindowListHead != NULL) {
    WindowListHead->wprev = newnode;
  }
  WindowListHead          = newnode;
  pop->hdc                = GetDC(pop->hWnd);
  
  pop->stringheight   = 10; 
  pop->stringwidth    = 6;
  pop->linewidth      = 1;   /* default pixel sizes of graphics until user changes them */
  pop->pointdiameter  = 1;
  pop->node           = newnode;
  
  newnode->bitwidth  = pop->w;
  newnode->bitheight = pop->h;
  
  /* Create and initialize primary graphics buffer */
  newnode->Buffer = CreateCompatibleDC(pop->hdc);
  newnode->BufferBit = CreateCompatibleBitmap(pop->hdc,pop->w,pop->h);
  newnode->store = SelectObject(newnode->Buffer,newnode->BufferBit);
  ExtFloodFill(newnode->Buffer,0,0,COLOR_WINDOW,FLOODFILLBORDER);
  
  
  newnode->event          = CreateEvent(NULL, TRUE, FALSE, NULL);
  newnode->DoubleBuffered = PETSC_FALSE;
  
  ReleaseDC(pop->hWnd,pop->hdc);
  ReleaseMutex(g_hWindowListMutex);
  PetscFunctionReturn(0);
}


EXTERN_C_BEGIN
int DrawCreate_Win32(Draw draw)
{       
  Draw_Win32  *windraw;
  HANDLE      hThread = NULL;
  int         ierr;
  WindowNode  newnode;
  
  windraw     = (Draw_Win32*)PetscMalloc(sizeof(Draw_Win32));CHKPTRQ(windraw);
  draw->data  = windraw;
  
  /* the following is temporary fix for initializing a global datastructure */
  if(!g_hWindowListMutex) {
    g_hWindowListMutex = CreateMutex(NULL,FALSE,NULL);
  }
  ierr = PetscMemcpy(draw->ops,&DvOps,sizeof(DvOps));CHKERRQ(ierr);
  
  windraw->hReadyEvent = CreateEvent(NULL,TRUE,FALSE,NULL);
  /* makes call to MessageLoopThread to creat window and attach a thread */
  CreateThread(NULL,0,(LPTHREAD_START_ROUTINE)MessageLoopThread_Win32,draw,0,(unsigned long*)hThread);
  CloseHandle(hThread);
  WaitForSingleObject(windraw->hReadyEvent,INFINITE);
  CloseHandle(windraw->hReadyEvent);
  WaitForSingleObject(g_hWindowListMutex,INFINITE);
  
  newnode                 = (WindowNode)PetscMalloc(sizeof(struct _p_WindowNode));CHKPTRQ(newnode);
  newnode->MouseListHead  = NULL;
  newnode->MouseListTail  = NULL;
  newnode->wnext          = WindowListHead;
  newnode->wprev          = NULL;
  newnode->hWnd           = windraw->hWnd;
  if(WindowListHead != NULL) {
    WindowListHead->wprev = newnode;
  }
  WindowListHead          = newnode;
  windraw->hdc            = GetDC(windraw->hWnd);
  
  windraw->stringheight   = 10; 
  windraw->stringwidth    = 6;
  windraw->linewidth      = 1;   /* default pixel sizes of graphics until user changes them */
  windraw->pointdiameter  = 1;
  windraw->node           = newnode;
  
  windraw->x = draw->x;
  windraw->y = draw->y;
  windraw->w = newnode->bitwidth    = draw->w;
  windraw->h = newnode->bitheight   = draw->h;  
  
  /* Create and initialize primary graphics buffer */
  newnode->Buffer = CreateCompatibleDC(windraw->hdc);
  newnode->BufferBit = CreateCompatibleBitmap(windraw->hdc,windraw->w,windraw->h);
  newnode->store = SelectObject(newnode->Buffer,newnode->BufferBit);
  ExtFloodFill(newnode->Buffer,0,0,COLOR_WINDOW,FLOODFILLBORDER);
  
  newnode->event          = CreateEvent(NULL,TRUE,FALSE,NULL);
  newnode->DoubleBuffered = PETSC_FALSE;
  
  ReleaseDC(windraw->hWnd,windraw->hdc);
  ReleaseMutex(g_hWindowListMutex);
  return 0;
}
EXTERN_C_END


/* FUNCTION: WndProc(HWND, unsigned, WORD, LONG)
   PURPOSE:  Processes messages for the main window.
   WM_COMMAND  - process the application menu
   WM_PAINT    - Paint the main window
   WM_DESTROY  - post a quit message and return */
LRESULT CALLBACK WndProc(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
{
  int         wmId, wmEvent;
  PAINTSTRUCT ps;
  
  switch (message) {
    HANDLE_MSG(hWnd,WM_PAINT,OnPaint_Win32);
    HANDLE_MSG(hWnd,WM_DESTROY,OnDestroy_Win32);
  case WM_COMMAND:
    wmId    = LOWORD(wParam); 
    wmEvent = HIWORD(wParam); 
    /* Parse the menu selections:*/
    switch (wmId) {
    case IDM_EXIT:
      DestroyWindow(hWnd);
      break;
    default:
      return DefWindowProc(hWnd, message, wParam, lParam);
    }
    break;
  case WM_LBUTTONUP:
    MouseRecord_Win32(hWnd,BUTTON_LEFT);
    break;
  case WM_RBUTTONUP:
    MouseRecord_Win32(hWnd,BUTTON_RIGHT);
    break;
  case WM_MBUTTONUP:
    MouseRecord_Win32(hWnd,BUTTON_CENTER);
    break; 
  default:
    return DefWindowProc(hWnd, message, wParam, lParam);
  }
  return 0;
}

static void OnPaint_Win32(HWND hWnd) 
{
  PAINTSTRUCT ps;
  HDC         hdc;
  WindowNode  current = NULL;
  InvalidateRect(hWnd,NULL,TRUE); 
  WaitForSingleObject(g_hWindowListMutex, INFINITE);
  current = WindowListHead;
  hdc     = BeginPaint(hWnd, &ps);
  
  while(current != NULL) {
    if (current->hWnd == hWnd) { 
      /* flushes primary buffer to window */
      BitBlt(hdc,
             0,0,
             GetDeviceCaps(hdc,HORZRES),
             GetDeviceCaps(hdc,VERTRES),
             current->Buffer,
             0,0,
             SRCCOPY);
      
      /* StretchBlt(hdc,
        0,0,
        w,h,
        current->Buffer,
        0,0,
        current->bitwidth,
        current->bitheight,
        SRCCOPY); */
      break;
    }
    current = current->wnext;
  }
  EndPaint(hWnd, &ps);
  ReleaseMutex(g_hWindowListMutex);
}

static int  MouseRecord_Win32(HWND hWnd,DrawButton button) 
{
        
  /* Called by all three mouse button actions
    Records needed mouse data in windows data structure */
  WindowNode current = NULL;
  MouseNode  newnode;
  POINT      mousepos;
  
  
  WaitForSingleObject(g_hWindowListMutex, INFINITE);
  current = WindowListHead;
  if(current->IsGetMouseOn == TRUE) {
    
    SetEvent(current->event);
    while (current != NULL) {   
      if(current->hWnd == hWnd) {       
        
        newnode = (MouseNode)PetscMalloc(sizeof(struct _p_MouseNode));CHKPTRQ(newnode);
        newnode->Button = button;
        GetCursorPos(&mousepos);
        newnode->user.x = mousepos.x;
        newnode->user.y = mousepos.y;
        ScreenToClient(hWnd,&mousepos);
        newnode->phys.x = mousepos.x;
        newnode->phys.y = mousepos.y;
        if(current->MouseListTail == NULL) {
          current->MouseListHead = newnode;
          current->MouseListTail = newnode;
        } else {
          current->MouseListTail->mnext = newnode;
          current->MouseListTail = newnode;
        }
        newnode->mnext = NULL;
        
        break;
      } 
      current = current->wnext;
    }
  }
  ReleaseMutex(g_hWindowListMutex);
  return 0;
}

static void OnDestroy_Win32(HWND hWnd) 
{
  /* searches linked list of window data and frees corresponding memory */
  WindowNode current;
  
  WaitForSingleObject(g_hWindowListMutex, INFINITE);
  current = WindowListHead;
  
  SetEvent(current->event);
  while (current != NULL) { 
    if(current->hWnd == hWnd) {
      if(current->wprev != NULL) {
        current->wprev->wnext = current->wnext;
      } else {
        WindowListHead = current->wnext;
      }
      if(current->MouseListHead) {
        deletemouselist_Win32(current);
      } else {
        PetscFree(current);
      }
      break;
    }
    current = current->wnext;
  }
  ReleaseMutex(g_hWindowListMutex);
  PostQuitMessage(0);
}
