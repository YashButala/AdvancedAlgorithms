#include<bits/stdc++.h>
#include "Voronoi_Gen.h"
using namespace std;
double getRadius(double x1,double y1,double	lx1, double ly1,double lx2,double ly2)
{
	double a=(ly1-ly2)/(lx2-lx1);
	double b=1;
	double c=-(lx2*a)-ly2;
	int flag=1;
	double radius= (a*x1 + b*y1 + c);
	double h=x1-(a*radius)/(a*a+b*b);
	double k=y1-(b*radius)/(a*a+b*b);
	if(h<=max(lx1,lx2)+0.0001 && h>=min(lx2,lx1)-0.0001)	
		if(k<=max(ly1,ly2)+0.0001 && k>=min(ly1,ly2)-0.0001)
		{	
			radius/=sqrt((a*a)+(b*b));
			flag=0;
		}
	 if(radius<0)
 			radius=-radius;
	if(flag)
		return min(min(abs(x1+250),abs(x1-250)),min(abs(y1+250),abs(y1-250)));
	else
		return min(radius,min(min(abs(x1+250),abs(x1-250)),min(abs(y1+250),abs(y1-250))));
	//	return radius;
}
void plot(Voronoi_Gen vdg,double x[],double y[], int n)
{
	double x1,y1,x2,y2;
    //fout in circles.svg file the x,y,co-ordinates of the center and radius of the circle.
    ofstream fout;
    fout.open("voronoi.svg");  
    int length=500;
    fout<<"<svg height=\""<<length<<"\" width=\""<<length<<"\">" <<endl;
    fout<<"<rect width=\""<<length<<"\" height=\""<<length<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />"<<endl;
    while(vdg.getNext(x1,y1,x2,y2))
	{
		printf("Line Segment of Voronoi (%f,%f)->(%f,%f)\n",x1,y1,x2, y2);
		fout<<"<line x1=\""<<x1+250<<"\" y1=\""<<y1+250<<"\" x2=\""<<x2+250<<"\" y2=\""<<y2+250<<"\" style=\"stroke:rgb(255,0,0);stroke-width:2\" />"<<endl;	
	}
	for(int i=0;i<n;i++)
	{
		fout<<"<circle cx=\""<<x[i]+250<<"\" cy=\""<<y[i]+250<<"\" r=\""<<2<<"\" stroke=\"black\" stroke-width=\"1\" fill=\"black\" fill-opacity=\"1.0\" />"<<endl;		
	}
 	for(int i=0;i<n;i++)
 	{
 		vdg.resetIterator();
 		vdg.getNext(x1,y1,x2,y2);
 		double r=getRadius(x[i],y[i],x1,y1,x2,y2);
 		while(vdg.getNext(x1,y1,x2,y2))
 			r=min(r,getRadius(x[i],y[i],x1,y1,x2,y2));
 		fout<<"<circle cx=\""<<x[i]+250<<"\" cy=\""<<y[i]+250<<"\" r=\""<<r<<"\" stroke=\"black\" stroke-width=\"1\" fill=\"blue\" fill-opacity=\"0.2\" />"<<endl;	
 	}
 	fout<<"</svg>"<<endl;
    fout.close(); 
}

int main(int argc,char **argv) 
{	
	srand(time(0));
	cout<<"Enter the number of points to be plotted"<<endl;
	int n;
	cin>>n;
	double xValues[n],yValues[n] ;
	for(int i=0;i<n;i++)
    {
    	bool flag=true;
        int a=(rand() % 500)-250;
        int b=(rand() % 500)-250;
        int c=(rand() % 500);
        int d=(rand() % 500);
        c=(a>=0?c:-c);
        d=(b>=0?d:-d);
        //cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
        double x=(double)a+(double)c/500.0;
        double y=(double)b+(double)d/500.0;
        for(int i1=0;i1<i;i1++)
        {
        	double dist=(x-xValues[i1])*(x-xValues[i1])+(y-yValues[i1])*(y-yValues[i1]);
        	if(dist<9.01)
        	{
        		i--;
        		flag=false;
        	}
        }
        if(flag)
        {	
       		xValues[i]=x;
        	yValues[i]=y;
        }	
    }
	
/*	xValues[0]=0;
	xValues[1]=30;
	xValues[2]=100;
	xValues[3]=2;
	xValues[4]=1;
	yValues[0]=0;
	yValues[1]=3;
	yValues[2]=100;
	yValues[3]=21;
	yValues[4]=30;*/
	long count = n;

	Voronoi_Gen vdg;
	vdg.generateVoronoi(xValues,yValues,count, -250,250,-250,250,3);

	vdg.resetIterator();

//	double x1,y1,x2,y2;

/*	printf("\n-------------------------------\n");
	while(vdg.getNext(x1,y1,x2,y2))
	{
		printf("GOT Line (%f,%f)->(%f,%f)\n",x1,y1,x2, y2);
		
	}
	vdg.resetIterator();*/
	plot(vdg,xValues,yValues,n);

	return 0;

	
}



