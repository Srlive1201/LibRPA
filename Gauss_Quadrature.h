#ifndef _GAUSS_QUADRATURE_H_
#define _GAUSS_QUADRATURE_H_

 #include <cstddef>


/*PURPOSE
   --Compute Gauss quadrature and its transformation from the unit interval [-1,1] to various type intervals and
     inverse transformation.

     Gauss Quadrature includes:
       1. Gauss-Legendre quadrature on [-1, 1].
       2. Gauss-Chebyshev quadrature of the first kind on [-1, 1].
       3. Gauss-Chebyshev quadrature of the second kind on [-1, 1].
 
     Transformation includes:
       1. [-1, 1]       -->   [x_one, x_two]
       2. [-1, 1]       -->   [0, inf)
       2. [-1, 1]       -->   (-inf, inf)
       3. [-1, 1]       -->   [x_zero, inf)
       4. [-1, 1]       -->   (-inf, x_zero)
       5. (-inf, inf)   -->   [-1, 1]
       6. [x_zero, inf)     -->   [-1, 1]
       7. (-inf, x_zero]    -->   [-1, 1]                     */
 class Gauss_Quadrature
 {

   public:

   char T{'L'};
   int N{40};
   double* nodes{NULL};
   double* weights{NULL};

/* Given a warning. Must be initialized EXPLICITLY! */
   Gauss_Quadrature( const char type = 'L', const int N = 40 );

   void reGenerate( const char type, const int N );

   void init( const char type, const int N );
   void free();

/* Query the functions of class Gauss_Quadrature. */
   void query();

/* Constructor.
     For T-type Gauss quadrature in the default/unit interval [-1,1]. Usually it is called by other types of constructors
     in the initialization of the T-type Gauss quadrature. */
   void unit( const char T = 'L' );

/* T-type Gauss quadrature in the interval [x_one,x_two]. */
   void close_Interval( const double x_one, const double x_two );

/* T-type Gauss quadrature in the interval [x_zero, infty] or [-infty,x_zero] which are specified by setting openside = 'R' and
   openside = 'L' respectively. */
   void open_Interval( const double x_zero, const char openside );

/* T-type Gauss quadrature in the interval (-infty,infty). */
   void infty_Interval( const char* inf_range );

/* Destructor. */
   virtual ~Gauss_Quadrature();


/* * Function Members * * * * * * * * * * * * * * * * * * * * * */

/*PURPOSE
   --Generate the Legendre polynomials and the corresponding derivatives.
     The weight function w(x) = 1. */
   static void Legendre_polynomial( const int N, const double x, double &y, double &dy );

/*PURPOSE
   *  --Generate the Gauss-Legendre Quadrature on the interval [-1,1].
   * N
   *  --The number of Gauss-Legendre grids.
   * GLq
   *  --Keep the Gauss-Legendre roots and corresponding weights (GLq[0][:],
   *    GLq[1][:]).
   */
   void Gaussq_Legendre(); 

/*PURPOSE
   --Generate the Gauss-Chebyshev Quadrature based on Chebyshev points of the first type (kind).
       The weght function w(x) = 1/ sqrt(1-x^2).
       The zeros of T_N(x) in [-1,1].
       ix = cos( (i+1/2)pi /N ), iwx = pi /N with, i = 0, 1, ..., N-1. */
   void Gaussq_ChebyshevT1();

/*PURPOSE
   --Get the Gauss-Chebyshev Quadrature based on Chebyshev points of the second type (kind), i.e., the nodes are the
     local extrema of T_N(x) in [-1,1].
       ix = cos( i*pi /(N+1) )
       iwx = pi/(N+1) *sin^2 ( i*pi /(N+1) ), with i = 1, ..., N. */
   void Gaussq_ChebyshevT2();

   void TrigontransGaussq();

/*PURPOSE
   --Linearly transforming Gauss Quadrature from [-1,1] to [x_zero,x_one].
        x' = (b-a)/2 *x + (a+b)/2
       wx' = (b-a)/2 *wx   */
   void transGaussq_linear( const double x_zero, const double x_one );

/*PURPOSE
   --Transforming the Gauss-Legendre Quadrature from [-1,1] to [0,inf). The modification proposed for the evaluation of
     the Casimir-Polder integral is applied.
       wx' = 1/(1-x)^2 *wx
        x' = 0.5* (1+x)/(1-x)             */
   void transGaussq_unit20_inf();

/*PURPOSE                                                               
   --Transforming Gauss Quadrature from [-1,1] to (-inf,inf).
       wx' = 5* (1+x^2)/(1-x^2)^2 *wx
        x' = 5* x/(1-x^2)    */
   void transGaussq_unit2inf();

/*PURPOSE
   --Transforming Gauss Quadrature from [-1,1] to [x_zero,inf).
       wx' = 30/(1-x)^2 *wx
        x' = 15* (x+1)/(1-x) +x_zero                         */
   void transGaussq_unit2x0_inf( const double x_zero );

/*PURPOSE
   --Transforming Gauss Quadrature from [-1,1] to (-inf,x_zero].
       wx' = 30 /(1+x)^2 *wx 
        x' = 15* (x-1)/(1+x) +x_zero                          */
   void transGaussq_unit2inf_x0( const double x_zero );
 
/*PURPOSE
   --Transforming the Gauss Quadrature from (-inf,inf) to [-1,1].
       wx' = wx* 10* ( 1+ 1/ sqrt( 25+ 4x^2 ) ) /( 5+ sqrt( 25+ 4x^2 ) )^2
        x' = 2x/ ( 5+ sqrt( 25+ 4x^2 ) )                                */
   void transGaussq_inf2unit();

/*PURPOSE
   --Transforming the Gauss Quadrature from [x_zero,inf) to [-1,1].
       wx' = wx* 30/(x+15-x_zero)^2
        x' = (x-15-x_zero)/(x+15-x_zero)                             */
   void transGaussq_x0_inf2unit( const double x_zero );

/*PURPOSE
   --Transforming the Gauss Quadrature from (-inf,x_zero] to [-1,1].
       wx' = wx* 30/(15-x+x_zero)^2
        x' = (15+x-x_zero)/(15-x+x_zero)                              */
   void transGaussq_inf_x02unit( const double x_zero );

 }; /* end of class Gauss_Quadrature. */


#endif 
/* _GAUSS_QUADRATURE_H_ */
