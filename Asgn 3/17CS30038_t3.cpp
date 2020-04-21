#include<bits/stdc++.h>
#include <fstream>
using namespace std;
struct point
{
	double x,y;
};
void mergex(point a[],int low ,int mid,int high)
{
    int i=low,j=mid+1;
    int k=0;
    point b[high-low+1];
    while(i<=mid&&j<=high)
    {
        if(a[i].x>a[j].x)
        {
            b[k]=a[j];
            j++;
            k++;
        }
        else
        {
            b[k]=a[i];
            i++;
            k++;
        }
    }
    while(i<=mid)
    {
        b[k]=a[i];
        i++;
        k++;
    }
    while(j<=high)
    {
        b[k]=a[j];
        j++;
        k++;
    }
    for(i=low,k=0;i<=high;i++,k++)
    {
        a[i]=b[k];
    }
    return ;
}
void xsort(point a[],int low,int high)
{
    if(low<high)
    {
        int mid=(low+high)/2;
        xsort(a,low,mid);
        xsort(a,mid+1,high);
        mergex(a,low,mid,high);
    }
}
void mergey(point a[],int low ,int mid,int high)
{
    int i=low,j=mid+1;
    int k=0;
    point b[high-low+1];
    while(i<=mid&&j<=high)
    {
        if(a[i].y>a[j].y)
        {
            b[k]=a[j];
            j++;
            k++;
        }
        else
        {
            b[k]=a[i];
            i++;
            k++;
        }
    }
    while(i<=mid)
    {
        b[k]=a[i];
        i++;
        k++;
    }
    while(j<=high)
    {
        b[k]=a[j];
        j++;
        k++;
    }
    for(i=low,k=0;i<=high;i++,k++)
    {
        a[i]=b[k];
    }
    return ;
}
void ysort(point a[],int low,int high)
{
    if(low<high)
    {
        int mid=(low+high)/2;
        ysort(a,low,mid);
        ysort(a,mid+1,high);
        mergey(a,low,mid,high);
    }
}
double dist(point a,point b)
{
	return sqrt((a.x-b.x)*(a.x-b.x)+(a.y-b.y)*(a.y-b.y));
}
double closestrlpair(point strip[],int n,double d)
{
	int i,j;
	double min=d;
  /*  for(int i=0;i<n;i++)
    {
        cout<<strip[i].x<<" "<<strip[i].y<<endl;
    }*/
	for(i=0;i<n;i++)
	{
		for(j=i+1;j<n && (strip[j].y-strip[i].y)<d;j++)
		{
			if(dist(strip[i],strip[j])<min)
			
				min=dist(strip[i],strip[j]);
	//	    cout<<i<<" "<<j<<" "<<dist(strip[i],strip[j])<<endl;  
		}
	}
	return min;
}
double closestpair(point x[],point y[],int l,int r)
{
	if (r-l==1)
	{
		return dist(x[l],x[r]);
	}
	else if(r-l==2)
	{
		double d1,d2,d3;
		d1=dist(x[l],x[l+1]);
		d2=dist(x[l+2],x[l+1]);
		d3=dist(x[l],x[l+2]);
		d1=(d1<d2?d1:d2);
		d1=(d1<d3?d1:d3);
		return d1;
				
	}
	else
	{
		int i,mid=(l+r)/2,yl=0,yr=0;
		double X=x[mid].x;
		point *y_l,*y_r;
		y_l=(point*)malloc((mid-l+1)*sizeof(point));
		y_r=(point*)malloc((r-mid)*sizeof(point));
        //divide using X coordinates
		for( i=0;i<=r-l;i++)
		{
			if(y[i].x<=X)
			{
				y_l[yl]=y[i];
				yl++;	
			}
			else
			{
				y_r[yr]=y[i];
				yr++;
			}
			
		}
		double d1=closestpair(x,y_l,l,mid);
		double d2=closestpair(x,y_r,mid+1,r);
		double d=(d1>d2?d2:d1);
   //     cout<<"d :"<<d<<endl;
		int n=0,j=0;
		for(i=0;i<=r-l;i++)
		{
			if(abs(y[i].x-X)<=d)
				n++;
		}
      //  cout<<"n :"<<n<<endl;
        //generate strip to find minimum in 2 seperations
		point *strip=(point*)malloc(sizeof(point)*n);
		for(i=0;i<=r-l;i++)
		{
			if(abs(y[i].x-X)<=d)
			{
				strip[j]=y[i];
				j++;
			}
		}		
		double d3=closestrlpair(strip,n,d);
		d3=(d3>d?d:d3);
		return d3;
	}
}
void generatecircles(point P[],int n,double k1)
{
    //fout in circles.svg file the x,y,co-ordinates of the center and radius of the circle.
    ofstream fout;
    fout.open("circles.svg");  
    int length=2*ceil(k1)+1000;
    fout<<"<svg height=\""<<length<<"\" width=\""<<length<<"\">" <<endl;
    fout<<"<rect width=\""<<length<<"\" height=\""<<length<<"\" style=\"fill:rgb(255,255,255); stroke-width:0; stroke:rgb(0,0,0)\" />"<<endl;
    for (int i = 0; i < n; ++i)
    { 
        fout << "<circle cx = \"" <<P[i].x +500 +k1<< "\" cy = \"" <<P[i].y +500+k1<< "\"  r = \"" << k1/2.0<<"\"" << " stroke=\"black\" stroke-width=\"1\" fill=\"blue\" fill-opacity=\"0.4\" />"<<endl;    
    } 
    fout<<"</svg>"<<endl;
    fout.close(); 
}
int main()
{ 
    int n;
    cout<<"Enter the number of points\n";
    cin>>n;
    if(n==1)
    {
        cout<<"Invalid!"<<endl;
        return 0;
    }
    point P[n];
    int sd;
    cout<<"Enter integer seed for random number generation\n";
    cin>>sd;
    srand(sd);


    //produce random numbers between -500.999 and +500.999 of floating type upto 3 decimals
    for(int i=0;i<n;i++)
    {
        int a=(rand() % 1000)-500;
        int b=(rand() % 1000)-500;
        int c=(rand() % 1000);
        int d=(rand() % 1000);
        c=(a>=0?c:-c);
        d=(b>=0?d:-d);
        //cout<<a<<" "<<b<<" "<<c<<" "<<d<<endl;
        float x=(float)a+(float)c/1000.0;
        float y=(float)b+(float)d/1000.0;
        P[i].x=x;
        P[i].y=y;
    }
    point Q[n];
    for(int i=0;i<n;i++)
        Q[i]=P[i];


//sort by x and y coordinates
    xsort(P,0,n-1);
    ysort(Q,0,n-1);

/*        for(int i=0;i<n;i++)
    {
        cout<<P[i].x<<" "<<P[i].y<<endl;
    }
*/
    double k1=closestpair(P,Q,0,n-1);
    cout<<"The smallest distance is "<<k1<<endl; 
    generatecircles(P,n,k1);
    return 0; 
 
}
