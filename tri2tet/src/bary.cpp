/*==============================================================================
 Programmed by: Christopher Simpson
 Date: 14 November 2014
 
 Details about Barycentric coordinates can be found on my blog at
 http://www.cdsimpson.net

 You are free to modify this for your purposes, and free to distribute.
 Please give me credit and link to my blog if you distribute.
*=============================================================================*/

void bary_tet( double *ans, double *p, double *a, double *b, double *c, double *d );
void bary_tri( double *ans, double *p, double *a, double *b, double *c );
void crossProduct(double *ans, double *v1, double *v2);
double dotProduct( double *v1, double *v2);

/*============================================================================*/
void bary_tet( double *ans, double *p, double *a, double *b, double *c, double *d )
/*==============================================================================
 This function gives the barycentric coordinates for a tetrahedron in 3 dimensions.
 It works by calculating the volume of the subtetrahedron using the scalar triple
 product. { V_tet = 1/6 * v1 * ( v2 x v2 ) }

 This works for points outside the tetrahedron as well.

 Input  Type[len]  Description
 -----  ---------  -----------
 *ans    double[4]  Pointer to an answer array. Values 0-3 will be overwritten
 *p      double[3]  Pointer to array of (x,y,z) coordinates for the test point
 *a      double[3]  Pointer to array of (x,y,z) coordinates for tet node 1
 *b      double[3]  Pointer to array of (x,y,z) coordinates for tet node 2
 *c      double[3]  Pointer to array of (x,y,z) coordinates for tet node 3
 *d      double[3]  Pointer to array of (x,y,z) coordinates for tet node 4

 No outputs other than *ans

 The node order for a tetrahedron is conterclockwise around the base, then up.

         4               3
          /|\             |.
         / | \            | . .
        /  |  \           |  .   .  2
       /   |   \          |    .  /\
      /    |    \         |     ./  \
     /     |     \        |     / .  \
    /      |      \       |    /   .  \
 0 /_______|_______\ 2    |   /     .  \
   \       |       /      |  /        . \
    \      |      /       | /          . \
     \     |     /        |/_____________.\
      \    |    /        0                  1
       \   |   /         
        \  |  /          
         \ | /           
          \|/            
            1            
         

*=============================================================================*/
{
  double vap[3];
  double vbp[3];
  double vcp[3];
  double vdp[3];
  double vab[3];
  double vac[3];
  double vad[3];
  double vbc[3];
  double vbd[3];
  double va;
  double vb;
  double vc;
  double vd;
  double v;
  double temp[3];

  int i;

  for ( i = 0; i < 3; i++ )
  {
    vap[i] = p[i] - a[i];
    vbp[i] = p[i] - b[i];
    vcp[i] = p[i] - c[i];
    vdp[i] = p[i] - d[i];
    vab[i] = b[i] - a[i];
    vac[i] = c[i] - a[i];
    vad[i] = d[i] - a[i];
    vbc[i] = c[i] - b[i];
    vbd[i] = d[i] - b[i];
  }
  crossProduct( temp, vbd, vbc );
  va = dotProduct( vbp, temp ) / 6.0 ;
  crossProduct( temp, vac, vad );
  vb = dotProduct( vap, temp ) / 6.0 ;
  crossProduct( temp, vad, vab );
  vc = dotProduct( vap, temp ) / 6.0 ;
  crossProduct( temp, vab, vac );
  vd = dotProduct( vap, temp ) / 6.0 ;
  crossProduct( temp, vac, vad );
  v  = dotProduct( vab, temp ) / 6.0 ;

  ans[0] = va / v;
  ans[1] = vb / v;
  ans[2] = vc / v;
  ans[3] = vd / v;
  
  return;
}

/*============================================================================*/
void bary_tri( double *ans, double *p, double *a, double *b, double *c )
/*==============================================================================
 This function gives the barycentric coordinates for a triangle in 3 dimensions.

 This works for points outside the triangle as well, and for triangles in any
 orientation. If the test point is not in the same plane as the triangle, the
 point is projected into the triangle plane, and the barycentric coordinates
 will point to the projected point.

 Input  Type[len]  Description
 -----  ---------  -----------
 *ans    double[3]  Pointer to an answer array. Values 0-2 will be overwritten
 *p      double[3]  Pointer to array of (x,y,z) coordinates for the test point
 *a      double[3]  Pointer to array of (x,y,z) coordinates for tet node 1
 *b      double[3]  Pointer to array of (x,y,z) coordinates for tet node 2
 *c      double[3]  Pointer to array of (x,y,z) coordinates for tet node 3

 No outputs other than *ans

 The nodes must be ordered counterclockwise around the triangle 
 (positive normal is v01 x v02)

        2
         /\
        /  \
       /    \
      /      \
     /        \
    /          \
   /            \
  /______________\
 0                1

*=============================================================================*/
{
  double vap[3];
  double vbp[3];
  double vcp[3];
  double vab[3];
  double vca[3];
  double vbc[3];
  double vac[3];
  double n[3];
  double na[3];
  double nb[3];
  double nc[3];

  int i;

  for ( i = 0; i < 3; i++ )
  {
    vap[i] = p[i] - a[i];
    vbp[i] = p[i] - b[i];
    vcp[i] = p[i] - c[i];

    vac[i] = c[i] - a[i];
    vab[i] = b[i] - a[i];
    vca[i] = a[i] - c[i];
    vbc[i] = c[i] - b[i];
  }

  crossProduct( n, vab, vac );
  crossProduct( na, vbc, vbp );
  crossProduct( nb, vca, vcp );
  crossProduct( nc, vab, vap );

  ans[0] = dotProduct( n, na ) / dotProduct( n, n );
  ans[1] = dotProduct( n, nb ) / dotProduct( n, n );
  ans[2] = dotProduct( n, nc ) / dotProduct( n, n );

  return;
}

/*============================================================================*/
void crossProduct(double *ans, double *v1, double *v2)
/*==============================================================================
 This function gives the cross product of two vectors.

 Input  Type[len]  Description
 -----  ---------  -----------
 *ans   double[3]  Pointer to an answer array. Values will be overwritten
 *v1    double[3]  Pointer to array of vector 1 components
 *v2    double[3]  Pointer to array of vector 2 components

 No outputs other than *ans

*=============================================================================*/
{
  ans[0] = v1[1]*v2[2] - v1[2]*v2[1];
  ans[1] = v1[2]*v2[0] - v1[0]*v2[2];
  ans[2] = v1[0]*v2[1] - v1[1]*v2[0];
  return;
}

/*============================================================================*/
double dotProduct( double *v1, double *v2 )
/*==============================================================================
 This function gives the dot product of two vectors.

 Input  Type[len]  Description
 -----  ---------  -----------
 *v1    double[3]  Pointer to array of vector 1 components
 *v2    double[3]  Pointer to array of vector 2 components

 Output  Type      Description
 ------ ------     -----------
 result double     the dot product

*=============================================================================*/
{
  double result = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
  return (result);
}

