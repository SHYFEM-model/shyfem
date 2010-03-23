
#include <stdio.h>
#include <string.h>
#include <grx.h>
#include <grdriver.h>

GR_DRIVER_MODE_ENTRY *ttable[];
GR_DRIVER_MODE_ENTRY *gtable[];

char buff[200];
char message[2000] = {""};

void main()

{
	int w,h;
	int col,white;
	int i,c=0;
	int colfree;
	int x,y,dx,dy,dxx,dyy;
	GR_DRIVER_MODE_ENTRY *tt;

/*	GrSetMode(GR_biggest_graphics); */
	GrSetMode(GR_default_graphics);
/*	GrSetMode(GR_80_25_text); */
/*	GrSetMode(GR_biggest_text); */
	GrGetDriverModes( ttable , gtable );


	printf("\nGraphic Modes :\n");
	sprintf(buff,"\nGraphic Modes :\n");
	strcat(message,buff);
	tt=*gtable;
	while(tt){
		w=tt->width;
		h=tt->height;
		c=tt->number_of_colors;
		if(w<=0) break;
		printf("%d %d %d\n",w,h,c);
		sprintf(buff,"%d %d %d\n",w,h,c);
		strcat(message,buff);
		tt++;
	}
	printf("\nText Modes :\n");
	sprintf(buff,"\nText Modes :\n");
	strcat(message,buff);
	tt=*ttable;
	while(tt){
		w=tt->width;
		h=tt->height;
		c=tt->number_of_colors;
		if(w<=0) break;
		printf("%d %d %d\n",w,h,c);
		sprintf(buff,"%d %d %d\n",w,h,c);
		strcat(message,buff);
		tt++;
	}


	w=GrScreenX();
	h=GrScreenY();
	col=GrNumColors();
	colfree=GrNumFreeColors();
	white=GrWhite();

	printf("\nScreen %d %d %d %d %d\n",w,h,col,colfree,white);
	sprintf(buff,"\nScreen %d %d %d %d %d\n",w,h,col,colfree,white);
	strcat(message,buff);

	c=getchar();

	GrClearScreen(5);
	GrLine(50,50,100,100,1);
	GrLine(150,150,100,100,white);

	dxx=dyy=25;
	x=dxx; y=200; dx=dy=10; 
	for(i=0;i<col;i++) {
	  if(x+dx+1>=w) {
		x=dxx;
		y+=dyy;
	  }
	  GrFilledBox(x-dx,y-dy,x+dx,y+dy,i);
	  x+=dxx;
	}

	c=getchar();
	
	GrSetMode(GR_80_25_text);
	
	puts(message);

}
