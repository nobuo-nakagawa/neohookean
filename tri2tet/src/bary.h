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

