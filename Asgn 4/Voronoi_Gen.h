#include <bits/stdc++.h>
#ifndef NULL
#define NULL 0
#endif
#define DELETED -2
#define le 0
#define re 1
struct	Freenode	
{
	struct	Freenode *nextfree;
};
typedef struct	Freenode FN;
struct FreeNodeArrayList
{
	struct	Freenode* memory;
	struct	FreeNodeArrayList* next;
};
typedef struct FreeNodeArrayList FNAL;
struct	Freelist	
{
	struct	Freenode	*head;
	int		nodesize;
};
typedef struct Freelist FL;
struct Point	
{
	double x,y;
};
typedef struct Point pt1;
struct Site 
{
	struct	Point	coord;
	int		sitenbr;
	int		refcnt;
};
typedef struct Site st;
struct Edge 	
{
	double   a,b,c;
	struct	Site 	*ep[2];
	struct	Site	*reg[2];
	int		edgenbr;

};
typedef struct Edge ed;
struct GraphEdge
{
	double x1,y1,x2,y2;
	struct GraphEdge* next;
};
typedef struct GraphEdge GEdge;
struct Halfedge
{
	struct	Halfedge	*Left_EL, *Right_EL;
	struct	Edge	*ELedge;
	int		ELrefcnt;
	char	ELpm;
	struct	Site	*vertex;
	double	ystar;
	struct	Halfedge *PQnext;
};
typedef struct Halfedge HEdge; 

int scomp(const void *p1,const void *p2)
{
	struct Point *s1 = (Point*)p1, *s2=(Point*)p2;
	if( (*s1).y < (*s2).y) return(-1);
	else if( (*s1).y > (*s2).y) return(1);
	else if((*s1).x < (*s2).x) return(-1);
	else if((*s1).x > (*s2).x) return(1);
	else return(0);
}
class Voronoi_Gen
{
public:
	Voronoi_Gen()
	{
		SiteIndex = 0;
		CellSite = 0;

		allMemoryList = new FNAL;
		(*allMemoryList).memory = 0;
		(*allMemoryList).next = 0;
		currentMemoryBlock = allMemoryList;
		allEdges = 0;
		iteratorEdges = 0;
		closestSiteDist = 0;
	}
	~Voronoi_Gen()
	{
		clear_all();
		remove_Edges();

		if(!allMemoryList) delete allMemoryList;
	}
	int generateVoronoi(double *x_coordinate, double *y_coordinate, int number_of_sites, double x_least, double x_max, double y_least, double y_max, double minimum_Dist=0)
	{
		clear_all();
		remove_Edges();


		closestSiteDist = minimum_Dist;

		NSites=number_of_sites; 
		init_erase(&sfl, sizeof (Site));
			
		CellSite = (st *) myalloc(NSites*sizeof( *CellSite));

		if(!CellSite)return 0;

		ymin = y_coordinate[0];
		xmax = x_coordinate[0];
		ymax = y_coordinate[0];
		x_min = x_coordinate[0];
		
		for(int j=0;j<NSites;j++)
		{
			CellSite[j].sitenbr = j;
			CellSite[j].refcnt = 0;	
			SiteIndex = 0;
		}
		int i = 0;
		while( i< NSites)
		{
			CellSite[i].coord.x = x_coordinate[i];
			CellSite[i].coord.y = y_coordinate[i];
			
			if(x_coordinate[i] < x_min)  	
				x_min = x_coordinate[i];
			else if(x_coordinate[i] > xmax) 
				xmax = x_coordinate[i];
			if(y_coordinate[i] < ymin) 	
				ymin = y_coordinate[i];
			else if(y_coordinate[i] > ymax)	
				ymax = y_coordinate[i];
			i++;
		}
		qsort(CellSite, NSites, sizeof (*CellSite), scomp);
		
		initialize_geo();
		double temp = 0;
		if(x_least > x_max)
		{
			x_least += x_max;
			x_max = x_least-x_max;
			x_least -= x_max;
		}
		if(y_least > y_max)
		{
			y_least += y_max;
			y_max = y_least-y_max;
			y_least -= y_max;
		}
		borderMinX = x_least;
		borderMaxX = x_max;
		borderMaxY = y_max;
		borderMinY = y_least;		
		SiteIndex = 0;
		triangulate = 0;
		voronoi(triangulate);

		return 1;
	}
	void resetIterator()
	{
		iteratorEdges = allEdges;
	}

	bool getNext(double& x1, double& y1, double& x2, double& y2)
	{
		if(iteratorEdges == 0)
			return false;
		
		else
		{ 
			x2 = (*iteratorEdges).x2;
			x1 = (*iteratorEdges).x1;
			y1 = (*iteratorEdges).y1;
			y2 = (*iteratorEdges).y2;
			iteratorEdges = iteratorEdges->next;
			return true;
		}	
	}


private:
	void clear_all()
	{
		if(CellSite != 0)
		{
			free(CellSite);
			CellSite = 0;
		}

		FNAL* current=0, *prev = 0;

		prev = allMemoryList;

		for(current = allMemoryList ;current->next != 0;)
		{
			prev = current;
			current = current->next;
			free(prev->memory);
			delete prev;
			prev = 0;
		}

		if(current != 0 && current->memory != 0)
		{
			free(current->memory);
			delete current;
		}

		allMemoryList = new FNAL;
		(*allMemoryList).next = 0;
		(*allMemoryList).memory = 0;
		currentMemoryBlock = allMemoryList;
	}
	void remove_Edges()
	{
		GraphEdge* geCurrent = 0, *gePrev = 0;
		geCurrent = gePrev = allEdges;

		while(geCurrent != 0 && (*geCurrent).next != 0)
		{
			gePrev = geCurrent;
			geCurrent = (*geCurrent).next;
			delete gePrev;
		}

		allEdges = 0;

	}
	char *getfree(FL *fl)
	{
		FN *t;

		if((*fl).head == (FN *) NULL)
		{	
			t =  (FN *) myalloc(sqrt_of_NSites * (*fl).nodesize);

			if(t == 0)
				return 0;
			
			(*currentMemoryBlock).next = new FNAL;
			currentMemoryBlock = currentMemoryBlock->next;
			(*currentMemoryBlock).memory = t;
			(*currentMemoryBlock).next = 0;

			for(int i=0; i<sqrt_of_NSites; i++) 	
				makefree((FN *)((char *)t+i*fl->nodesize), fl);		
		};
		t = (*fl).head;
		(*fl).head = ((*fl).head) -> nextfree;
		return((char *)t);
	}

	HEdge *PQfind();
	int PQempty()
	{
		return(seg_count==0);
	}
	
	HEdge **Hash_EL;
	HEdge *HEcreate(), *Left_EL();
	HEdge *Right_EL(), *ELleftbnd();
	HEdge *HEcreate(ed  *e,int pm)
	{
		HEdge *answer;
		answer = (HEdge *) getfree(&hfl);
		(*answer).ELedge = e;
		(*answer).ELpm = pm;
		(*answer).PQnext = (HEdge *) NULL;
		(*answer).vertex = (st *) NULL;
		(*answer).ELrefcnt = 0;
		return(answer);
	}


	struct Point PQ_min()
	{
		struct Point answer;
		while(hash_Seg[seg_min].PQnext == (HEdge *)NULL) 
		{
			seg_min ++;
		};
		answer.x = (*(hash_Seg[seg_min].PQnext)).vertex -> coord.x;
		answer.y = (*(hash_Seg[seg_min].PQnext)).ystar;
		return (answer);
	}
	HEdge *PQextractmin()
	{
		HEdge *curr;
		
		curr = hash_Seg[seg_min].PQnext;
		hash_Seg[seg_min].PQnext = (*curr).PQnext;
		seg_count --;return(curr);
	}
	void init_erase(FL *fl,int size)
	{
		(*fl).head = (FN *) NULL;
		(*fl).nodesize = size;
	}
	void makefree(FN *curr,FL *fl)
	{
		(*curr).nextfree = (*fl).head;
		(*fl).head = curr;
	}
	void initialize_geo()
	{	
		double sn;
		init_erase(&efl, sizeof(Edge));
		No_of_vertices = 0;
		nedges = 0;
		sn = (double)NSites+4;
		sqrt_of_NSites = (int)sqrt(sn);
		deltay = ymax;
		deltay-= ymin;
		deltax = xmax;
		deltay-= x_min;
	}
	bool voronoi(int triangulate)
	{
		st *newsite, *down, *up, *temp, *p;
		st *v;
		struct Point newintstar;
		int pm;
		HEdge *lbnd, *rbnd, *llbnd, *rrbnd, *bisector;
		ed  *e;
		
		PQinitialize();
		bottomsite = nextone();
		bool retval = ELinitialize();

		if(!retval) return false;
		
		newsite = nextone();
		while(1)
		{

			if(!PQempty()) newintstar = PQ_min();
			
			//if the lowest site has a smaller y value than the lowest vector intersection, process the site
			//otherwise process the vector intersection		

			if (newsite != (st *)NULL 	&& (PQempty() || newsite -> coord.y < newintstar.y
				|| (newsite->coord.y == newintstar.y && newsite->coord.x < newintstar.x)))
			{				
				lbnd = ELleftbnd(&(newsite->coord));				//get the first HalfEdge to the LEFT of the new site
				rbnd = Right_EL(lbnd);						//get the first HalfEdge to the RIGHT of the new site
				down = rightreg(lbnd);						//if this halfedge has no edge, , down = bottom site (whatever that is)
				e = bisect(down, newsite);					//create a new edge that bisects 
				bisector = HEcreate(e, le);					//create a new HalfEdge, setting its ELpm field to 0			
				ELinsert(lbnd, bisector);					//insert this new bisector edge between the left and right vectors in a linked list	

				if ((p = intersect(lbnd, bisector)) != (st *) NULL) 	//if the new bisector intersects with the left edge, 
				{	
					seg_delete(lbnd);							//remove the left edge's vertex, and put in the new one
					seg_insert(lbnd, p, dist(p,newsite));
				};
				lbnd = bisector;						
				bisector = HEcreate(e, re);					//create a new HalfEdge, setting its ELpm field to 1
				ELinsert(lbnd, bisector);					//insert the new HE to the right of the original bisector earlier in the IF stmt

				if ((p = intersect(bisector, rbnd)) != (st *) NULL)	//if this new bisector intersects with the
				{	
					seg_insert(bisector, p, dist(p,newsite));			//push the HE into the ordered linked list of vertices
				};
				newsite = nextone();	
			}
			else if (!PQempty()) /* intersection is smallest - this is a vector event */			
			{	
				lbnd = PQextractmin();						//pop the HalfEdge with the lowest vector off the ordered list of vectors				
				llbnd = Left_EL(lbnd);						//get the HalfEdge to the left of the head HE
				rbnd = Right_EL(lbnd);						//get the HalfEdge to the right of the head HE
				rrbnd = Right_EL(rbnd);						//get the HalfEdge to the right of the HE to the right of the lowest HE 
				down = leftreg(lbnd);						//get the Site to the left of the left HE which it bisects
				up = rightreg(rbnd);						//get the Site to the right of the right HE which it bisects

			
				v = (*lbnd).vertex;						//get the vertex that caused this event
				makevertex(v);							//set the vertex number - couldn't do this earlier since we didn't know when it would be processed
				set_vertex_end((*lbnd).ELedge,(*lbnd).ELpm,v);	//set the endpoint of the left HalfEdge to be this vector
				set_vertex_end((*rbnd).ELedge,(*rbnd).ELpm,v);	//set the endpoint of the right HalfEdge to be this vector
				ELdelete(lbnd);							//mark the lowest HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
				seg_delete(rbnd);							//remove all vertex events to do with the  right HE
				ELdelete(rbnd);							//mark the right HE for deletion - can't delete yet because there might be pointers to it in Hash Map	
				pm = le;								//set the pm variable to zero
				
				if ((*down).coord.y > (*up).coord.y)		//if the site to the left of the event is higher than the Site
				{										//to the right of it, then swap them and set the 'pm' variable to 1
					temp = down; 
					down = up; 
					up = temp; 
					pm = re;
				}
				e = bisect(down, up);					//create an Edge (or line) that is between the two CellSite. This creates
														//the formula of the line, and assigns a line number to it
				bisector = HEcreate(e, pm);				//create a HE from the Edge 'e', and make it point to that edge with its ELedge field
				ELinsert(llbnd, bisector);				//insert the new bisector to the right of the left HE
				set_vertex_end(e, re-pm, v);					//set one set_vertex_end to the new edge to be the vector point 'v'.
														//If the site to the left of this bisector is higher than the right
														//Site, then this set_vertex_end is put in position 0; otherwise in pos 1
				del_vec(v);								//delete the vector 'v'

				//if left HE and the new bisector don't intersect, then delete the left HE, and reinsert it 
				if((p = intersect(llbnd, bisector)) != (st *) NULL)
				{	
					seg_delete(llbnd);
					seg_insert(llbnd, p, dist(p,down));
				};

				//if right HE and the new bisector don't intersect, then reinsert it 
				if ((p = intersect(bisector, rrbnd)) != (st *) NULL)
				{	
					seg_insert(bisector, p, dist(p,down));
				};
			}
			else 
				break;
		};

		


		for(lbnd=Right_EL(left_end_EL); lbnd != right_end_EL; lbnd=Right_EL(lbnd))
		{	
			e = (*lbnd).ELedge;
			clip_line(e);
		};

		clear_all();

		return true;
		
	}
	void ref(st *v)
	{
		v -> refcnt += 1;
	}
	void del_vec(st *v)
	{
		v -> refcnt -= 1;
		if (v -> refcnt == 0 ) 
			makefree((Freenode*)v, &sfl);
	}
	void set_vertex_end(ed  *e,int lr,st * s)
	{
		e -> ep[lr] = s;
		ref(s);
		if(e -> ep[re-lr]== (st *) NULL) 
			return;

		clip_line(e);

		del_vec((*e).reg[le]);
		del_vec((*e).reg[re]);
		makefree((Freenode*)e, &efl);
	}
	void ELdelete(HEdge *he)
	{
		(he -> Left_EL) -> Right_EL = he -> Right_EL;
		(he -> Right_EL) -> Left_EL = he -> Left_EL;
		he -> ELedge = (ed  *)DELETED;
	}	

	HEdge *ELleftbnd(struct Point *p)
	{
		int i, bucket;
		HEdge *he;
		
		/* Use hash table to get close to desired halfedge */
		bucket = (int)(((*p).x - x_min)/deltax * ELhashsize);	//use the hash function to find the place in the hash map that this HalfEdge should be

		if(bucket<0) bucket =0;					//make sure that the bucket position in within the range of the hash array
		if(bucket>=ELhashsize) bucket = ELhashsize - 1;

		he = ELgethash(bucket);
		if(he == (HEdge *) NULL)			//if the HE isn't found, search backwards and forwards in the hash map for the first non-null entry
		{   
			for(i=1; 1 ; i += 1)
			{	
				if ((he=ELgethash(bucket-i)) != (HEdge *) NULL) 
					break;
				if ((he=ELgethash(bucket+i)) != (HEdge *) NULL) 
					break;
			};
		};
	//	ntry += 1;
		/* Now search linear list of halfedges for the correct one */
		if (he==left_end_EL  || (he != right_end_EL && right_of(he,p)))
		{
			do 
			{
				he = (*he).Right_EL;
			} 
			while (he!=right_end_EL && right_of(he,p));	//keep going right on the list until either the end is reached, or you find the 1st edge which the point
			he = (*he).Left_EL;				//isn't to the right of
		}
		else 							//if the point is to the left of the HalfEdge, then search left for the HE just to the left of the point
			do 
			{
				he = he -> Left_EL;
			} 
			while (he!=left_end_EL && !right_of(he,p));
			
		/* Update hash table and reference counts */
		if(bucket > 0 && bucket+1 <ELhashsize)
		{	
			if(Hash_EL[bucket] != (HEdge *) NULL) 
			{
				(*Hash_EL[bucket] ).ELrefcnt -= 1;
			}
			Hash_EL[bucket] = he;
			(*Hash_EL[bucket]).ELrefcnt += 1;
		}
		return (he);
	}
	HEdge *Right_EL(HEdge *he)
	{
		return ((*he).Right_EL);
	}
	void makevertex(st *v)
	{
		(*v).sitenbr = No_of_vertices;
		No_of_vertices++;
	}

	void seg_insert(HEdge *he,st * v, double offset)
	{
		HEdge *last, *next;
		
		he -> vertex = v;
		ref(v);
		he -> ystar = (double)(v -> coord.y + offset);
		last = &hash_Seg[PQbucket(he)];
		while ((next = last -> PQnext) != (HEdge *) NULL &&
			((*he).ystar  > (*next).ystar  ||
			((*he).ystar == (*next).ystar && (*v).coord.x > (*next).vertex->coord.x)))
		{	
			last = next;
		};
		(*he).PQnext = (*last).PQnext; 
		(*last).PQnext = he;
		seg_count +=1;
	}
	void seg_delete(HEdge *he)
	{
		HEdge *last;
		
		if( (*he).vertex != (st *) NULL)
		{	
			last = &hash_Seg[PQbucket(he)];
			while ((*last).PQnext != he) 
				last = (*last).PQnext;

			(*last).PQnext = (*he).PQnext;
			seg_count--;
			del_vec(he -> vertex);
			he -> vertex = (st *) NULL;
		};
	}


	bool ELinitialize()
	{
		int i;
		init_erase(&hfl, sizeof **Hash_EL);
		ELhashsize = 2 * sqrt_of_NSites;
		Hash_EL = (HEdge **) myalloc ( sizeof *Hash_EL * ELhashsize);

		if(Hash_EL == 0)
			return false;

		for(i=0; i<ELhashsize; i++) 
			Hash_EL[i] = (HEdge *)NULL;
		left_end_EL = HEcreate( (ed  *)NULL, 0);
		right_end_EL = HEcreate( (ed  *)NULL, 0);
		left_end_EL -> Left_EL = (HEdge *)NULL;
		left_end_EL -> Right_EL = right_end_EL;
		right_end_EL -> Left_EL = left_end_EL;
		right_end_EL -> Right_EL = (HEdge *)NULL;
		Hash_EL[0] = left_end_EL;
		Hash_EL[ELhashsize-1] = right_end_EL;

		return true;
	}
	void ELinsert(struct	Halfedge *lb, HEdge *newHe)
	{
		(*newHe).Left_EL = lb;
		(*newHe).Right_EL = lb -> Right_EL;
		(lb -> Right_EL) -> Left_EL = newHe;
		(*lb).Right_EL = newHe;
	}
	HEdge * ELgethash(int b)
	{
		HEdge *he;
		
		if(b<0 || b>=ELhashsize) 
			return((HEdge *) NULL);
		he = Hash_EL[b]; 
		if (he == (HEdge *) NULL || he->ELedge != (ed  *) DELETED ) 
			return (he);
		
		/* Hash table points to deleted half edge.  Patch as necessary. */
		Hash_EL[b] = (HEdge *) NULL;
		if ((he -> ELrefcnt -= 1) == 0) 
			makefree((Freenode*)he, &hfl);
		return ((HEdge *) NULL);
	}	
	HEdge *Left_EL(HEdge *he)
	{
		return (he -> Left_EL);
	}
	st *leftreg(HEdge *he)
	{
		if(he -> ELedge == (ed  *)NULL) 
			return(bottomsite);
		return( he -> ELpm == le ? 
			he -> ELedge -> reg[le] : he -> ELedge -> reg[re]);
	}

	bool PQinitialize()
	{
		int i; 
		
		seg_count = 0;
		seg_min = 0;
		hash_seg_size = 4 * sqrt_of_NSites;
		hash_Seg = (HEdge *) myalloc(hash_seg_size * sizeof *hash_Seg);

		if(hash_Seg == 0)
			return false;

		for(i=0; i<hash_seg_size; i+=1) hash_Seg[i].PQnext = (HEdge *)NULL;

		return true;
	}
	int PQbucket(HEdge *he)
	{
		int bucket;
		
		bucket = (int)((he->ystar - ymin)/deltay * hash_seg_size);
		if (bucket<0) bucket = 0;
		if (bucket>=hash_seg_size) bucket = hash_seg_size-1 ;
		if (bucket < seg_min) seg_min = bucket;
		return(bucket);
	}
	void clip_line(ed  *e)
	{
		st *s1;
		double x1=0,x2=0,y1=0,y2=0, temp = 0;;

		x1 = (*((*e).reg[0])).coord.x;
		x2 = (*((*e).reg[1])).coord.x;
		y1 = (*((*e).reg[0])).coord.y;
		y2 = (*((*e).reg[1])).coord.y;
		double k=(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1)));
		k=sqrt(k);
		//if the distance between the two points this line was created from is less than 
		//the that between closest CellSite, then ignore it
		if( k< closestSiteDist)
		{
			return;
		}
		xmin_pt = borderMinX;x_max_pt = borderMaxX;
		y_min_pt = borderMinY;y_max_pt = borderMaxY;
		st *s2;
		if(e -> a == 1.0 && e ->b >= 0.0)
		{	
			s1 = (*e).ep[1];
			s2 = (*e).ep[0];
		}
		else 
		{
			s1 = (*e).ep[0];
			s2 = (*e).ep[1];
		};
		
		if((*e).a == 1.0)
		{
			y1 = y_min_pt;
			if (s1!=(st *)NULL && s1->coord.y > y_min_pt)
			{
				y1 = s1->coord.y;
			}
			if(y1>y_max_pt) 
			{
			//	printf("\nClipped (1) y1 = %f to %f",y1,y_max_pt);
				y1 = y_max_pt;
				//return;
			}
			x1 = (*e).c - (*e).b * y1;
			y2 = y_max_pt;
			if (s2!=(st *)NULL && (*s2).coord.y < y_max_pt) 
				y2 = (*s2).coord.y;

			if(y2<y_min_pt) 
			{
				//printf("\nClipped (2) y2 = %f to %f",y2,y_min_pt);
				y2 = y_min_pt;
				//return;
			}
			x2 = ((*e).c) - ((*e).b) * y2;
			if (((x1> x_max_pt) & (x2>x_max_pt)) | ((x1<xmin_pt)&(x2<xmin_pt))) 
			{
				//printf("\nClipLine jumping out(3), x1 = %f, xmin_pt = %f, x_max_pt = %f",x1,xmin_pt,x_max_pt);
				return;
			}
			if(x1> x_max_pt)
			{	
				x1 = x_max_pt; 
				y1 = ((*e).c - x1)/((*e).b);
			}
			if(x1<xmin_pt)
			{
				x1 = xmin_pt; 
				y1 = ((*e).c - x1)/((*e).b);
			}
			if(x2>x_max_pt)
			{	
				x2 = x_max_pt; 
				y2 = ((*e).c - x2)/((*e).b);
			}
			if(x2<xmin_pt)
			{	
				x2 = xmin_pt;
				y2 = ((*e).c - x2)/((*e).b);
			}
		}
		else
		{
			x1 = xmin_pt;
			if (s1!=(st *)NULL && s1->coord.x > xmin_pt) 
				x1 = s1->coord.x;
			if(x1>x_max_pt) 
			{
				//printf("\nClipped (3) x1 = %f to %f",x1,xmin_pt);
				//return;
				x1 = x_max_pt;
			}
			y1 = e -> c - e -> a * x1;
			x2 = x_max_pt;
			if (s2!=(st *)NULL && s2->coord.x < x_max_pt) 
				x2 = (*s2).coord.x;
			if(x2<xmin_pt) 
			{
				//printf("\nClipped (4) x2 = %f to %f",x2,xmin_pt);
				//return;
				x2 = xmin_pt;
			}
			y2 = (*e).c - (*e).a * x2;
			if (((y1> y_max_pt) & (y2>y_max_pt)) | ((y1<y_min_pt)&(y2<y_min_pt))) 
			{
				//printf("\nClipLine jumping out(6), y1 = %f, y_min_pt = %f, y_max_pt = %f",y2,y_min_pt,y_max_pt);
				return;
			}
			if(y1> y_max_pt)
			{	
				y1 = y_max_pt;
				x1 = ((*e).c - y1)/(*e).a;
			}
			if(y1<y_min_pt)
			{	
				y1 = y_min_pt; 
				x1 = ((*e).c - y1)/(*e).a;
			}
			if(y2>y_max_pt)
			{	y2 = y_max_pt; 
				x2 = ((*e).c - y2)/(*e).a;
			}
			if(y2<y_min_pt)
			{	
				y2 = y_min_pt; 
				x2 = ((*e).c - y2)/(*e).a;
			}
		}
		
		//printf("\nPushing line (%f,%f,%f,%f)",x1,y1,x2,y2);
		line(x1,y1,x2,y2);
	}

	char *myalloc(unsigned n)
	{
		char *t=(char*)malloc(n);
//		total_alloc += n;
		return(t);
	}
	int right_of(HEdge *el,struct Point *p)
	{
		ed  *e;
		st *topsite;
		int right_of_site, head, fast;
		double dxp, dyp;
		e = (*el).ELedge;
		topsite = (*e).reg[1];
		if( (*p).x > (*topsite).coord.x)
			right_of_site=true;
		else
			right_of_site=false;
		if(right_of_site && (*el).ELpm == le) 
			return 1;
		if(!right_of_site && (*el).ELpm == re) 
			return 0;
		double  t1, t2,t3;
		double  dxs, yl;
		if ((*e).a == 1.0)
		{	
			dyp = (*p).y; 

			dxp = (*p).x ;
			fast = 0;
			dxp -= (*topsite).coord.x;
			dyp-=topsite->coord.y;
			if ((!right_of_site & ((*e).b<0.0)) | (right_of_site & ((*e).b>=0.0)) )
			{	
				head = dyp>= (*e).b*dxp;	
				fast = head;
			}
			else 
			{	
				head = p->x + p->y*(*e).b > (*e).c;
				if((*e).b<0.0) 
					head = !head;
				if (!head) 
					fast = 1;
			};
			if (!fast)
			{	
				dxs = (*topsite).coord.x ;
				dxs-= ((*e).reg[0])->coord.x;
				head = (*e).b * (dxp*dxp - dyp*dyp) < dxs*dyp*(1.0+2.0*dxp/dxs + (*e).b*(*e).b);
				if((*e).b<0.0) 
					head = !head;
			}
		}
		else  /*e->b==1.0 */
		{	
			yl = (*e).c; 
			yl -= (*e).a*p->x;
			t1 = (*p).y - yl;
			t2 = (*p).x - (*topsite).coord.x;
			t3 = yl - (*topsite).coord.y;
			head = t1*t1 > t2*t2 + t3*t3;
		}
		return ((*el).ELpm==le ? head : !head);
	}

	st *rightreg(HEdge *he)
	{
		if(he -> ELedge == (ed  *)NULL) //if this halfedge has no edge, return the bottom site (whatever that is)
			return(bottomsite);

		//if the ELpm field is zero, return the site 0 that this edge bisects, otherwise return site number 1
		return( he -> ELpm == le ? he -> ELedge -> reg[re] : he -> ELedge -> reg[le]);
	}
	ed  *bisect(struct	Site *s1,struct	Site *s2)
	{
		double dx,dy;

		ed  *newedge;	

		newedge = (ed  *) getfree(&efl);
		
		newedge -> reg[0] = s1; //store the CellSite that this edge is bisecting
		newedge -> reg[1] = s2;
		ref(s1); 
		ref(s2);
		newedge -> ep[0] = (st *) NULL; //to begin with, there are no endpoints on the bisector - it goes to infinity
		newedge -> ep[1] = (st *) NULL;
		
		dx = s2->coord.x - s1->coord.x;			//get the difference in x dist between the CellSite
		dy = s2->coord.y - s1->coord.y;
		double adx,ady;
		if(dx>0)					//make sure that the difference in positive
			adx=dx;
		else
			adx=-dx;
		if(dy>0)					//make sure that the difference in positive
			ady=dy;
		else
			ady=-dy;						
		newedge -> c = (double)(s1->coord.x * dx + s1->coord.y * dy + (dx*dx + dy*dy)*0.5);//get the slope of the line

		if (adx>ady)
		{	
			(*newedge).a = 1.0; 
			(*newedge).b = dy/dx; 
			(*newedge).c /= dx;//set formula of line, with x fixed to 1
		}
		else
		{	
			(*newedge).b = 1.0; 
			(*newedge).a = dx/dy; 
			(*newedge).c /= dy;//set formula of line, with y fixed to 1
		};
		
		(*newedge).edgenbr = nedges;

		//printf("\nbisect(%d) ((%f,%f) and (%f,%f)",nedges,s1->coord.x,s1->coord.y,s2->coord.x,s2->coord.y);
		
		nedges ++;
		return(newedge);
	}
	double dist(st *s,st *t)
	{
		double dx,dy;
		dx = (*s).coord.x;
		dx -= (*t).coord.x;
		dy = (*s).coord.y ;
		dy-= (*t).coord.y;
		return (double)(sqrt(dx*dx + dy*dy));
	}

	st *intersect(HEdge *el1, HEdge *el2, struct Point *p=0)
	{
		struct	Edge *e1,*e2, *e;
		struct  Halfedge *el;
		double d, xint, yint;
		int right_of_site;
		st *v;
		
		e1 = el1 -> ELedge;
		e2 = el2 -> ELedge;
		if(e1 == (ed *)NULL || e2 == (ed *)NULL) 
			return ((st *) NULL);

		//if the two edges bisect the same parent, return null
		if ((*e1).reg[1] == (*e2).reg[1]) 
			return ((st *) NULL);
		
		d = e1->a * e2->b - e1->b * e2->a;
		if (-1.0e-10<d && d<1.0e-10) 
			return ((st *) NULL);
		
		xint = (e1->c*e2->b - e2->c*e1->b);
		yint = (e2->c*e1->a - e1->c*e2->a);
		xint /= d;
		yint /=d;
		
		if( ((*e1).reg[1]->coord.y < (*e2).reg[1]->coord.y) ||((*e1).reg[1]->coord.y == (*e2).reg[1]->coord.y && (*e1).reg[1]->coord.x < (*e2).reg[1]->coord.x) )
		{	
			el = el1; 
			e = e1;
		}
		else
		{	
			el = el2; 
			e = e2;
		};
		
		right_of_site = xint >= (*((*e).reg[1])).coord.x;
		if ((right_of_site && el -> ELpm == le) || (!right_of_site && (*el).ELpm == re)) 
			return ((st *) NULL);
		
		//create a new site at the point of intersection - this is a new vector event waiting to happen
		v = (st *) getfree(&sfl);
		 (*v).refcnt = 0;
		(*v).coord.x = xint;
		(*v).coord.y = yint;
		return(v);
	}


	st *nextone()
	{
		st *s;
		if(SiteIndex < NSites)
		{	
			s = &CellSite[SiteIndex];
			SiteIndex += 1;
			return(s);
		}
		else
		{
			return( (st *)NULL);
		}	
			
	}
	void pushGraphEdge(double x1, double y1, double x2, double y2)
	{
		GraphEdge* newEdge = new GraphEdge;
		(*newEdge).next = allEdges;
		allEdges = newEdge;
		(*newEdge).x1 = x1;
		(*newEdge).y1 = y1;
		(*newEdge).x2 = x2;
		(*newEdge).y2 = y2;
	}

//	void openpl()

	void line(double pt1_x, double pt1_y, double pt2_x, double pt2_y)
	{	
		pushGraphEdge(pt1_x,pt1_y,pt2_x,pt2_y);

	}
//	void circle(double x, double y, double radius){}

	FL	hfl;
	HEdge *left_end_EL, *right_end_EL;
	int 	ELhashsize, triangulate, NSites, SiteIndex, sqrt_of_NSites, No_of_vertices; 
	double	x_min, xmax, ymin, ymax, closestSiteDist;
	struct	Site	*CellSite;
	FL sfl;
	struct	Site	*bottomsite;

	int		nedges;
	FL efl;
	int		hash_seg_size;
	HEdge *hash_Seg;
	int		seg_count;
	int		seg_min;

		

	double xmin_pt, x_max_pt, y_min_pt, y_max_pt, borderMinX, borderMaxX, borderMinY, borderMaxY;

	FNAL* allMemoryList;
	FNAL* currentMemoryBlock;

	GraphEdge* allEdges;
	GraphEdge* iteratorEdges;

	double deltax, deltay;
	
};





