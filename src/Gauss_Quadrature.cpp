 #include "Gauss_Quadrature.h"
 //#include "../src_global/constants.h"
 #include "input.h"
 #include "constants.h"
 #include <cmath>
 #include <iostream>
 #include <iomanip>
 #include <cstddef>
using namespace std;
/*Constructor.
   --Gives a warning when initialized implicitly. */
 Gauss_Quadrature::Gauss_Quadrature( const char type, const int n_points )
 {
   init( type, n_points );
 }

 void Gauss_Quadrature::reGenerate( const char type, const int n_points )
 {
   free();
   init( type, n_points );

   return;
 }

/*Constructor.
   --Gives a warning when initialized implicitly. */
 void Gauss_Quadrature::init( const char type, const int n_points )
 {
   T = type;
   N = n_points;
   nodes = new double [ N ]{0};
   weights = new double [ N ]{0};
   std::cout << "  | Gauss quadrature type:  " << T << "\n    The number of points: " << std::setw(3) << N << '\n';

   return;
 }

/*Constructor.
   --Query the functions of class Gauss_Quadrature. */
 void Gauss_Quadrature::query()
 {
   std::cout << "  | Class \"Gauss_quadrature\" is designed for generating Gauss quadrature.\n\n";
   std::cout << "    Gauss-Legendre quadrature, Gauss-Chebyshev quadrature of the first and the second kind are " \
                "available.\n Range choices are [x_one, x_two], [x_zero, infty), (-infty,x_zero] and (-infty,infty).\n";
   return;
 }

/*Constructor.
   --For T-type Gauss quadrature in the default/initial interval [-1,1]. */
 void Gauss_Quadrature::unit( const char T )
 {
   switch(T)
   {
     case 'L':
       std::cout << "    Generates Gauss-Legendre quadrature ";
       Gaussq_Legendre();
       break;
     case 'C':
       std::cout << "    Generates Gauss-Chebyshev quadrature of the first kind ";
       Gaussq_ChebyshevT1();
       break;   
     case 'U':
       std::cout << "    Generates Gauss-Chebyshev quadrature of the second kind ";
       Gaussq_ChebyshevT2();
       break;
     default:
       std::cout << "\n   ! Wrong assignment for the Gauss Quadrature type!\n\n Available choices:\n\t'L\': " \
                    "Gauss-Legendre quadrature.\n\t'C\': Gauss-Chebyshev quadrature of the first kind."       \
                    "\n\t\'U\': Gauss-Chebyshev quadrature of the second kind ";
       break;
   }
   std::cout << "in the interval: [-1.00, 1.00]\n";

   return;
 }

/*Constructor.
   --Generates T-type Gauss quadrature in the interval [x_one,x_two]. */
 void Gauss_Quadrature::close_Interval( const double x_one, const double x_two )
 {
   if( x_two <= x_one )
   {
     std::cout << "  ! ERROR: Sick range! Check the specified interval.\n";
     exit(1);
   }
   unit(T);
   std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(1)\
             << "    Transforms to [" << std::setprecision(2) << x_one << ", " << x_two << "].\n\n";
   transGaussq_linear( x_one, x_two );

   return;
 }

/*Contructor.
   --Generates T-type Gauss quadrature in the interval [x_zero, infty] or [-infty,x_zero] which are specified by setting
     openside = 'R' and openside = 'L' respectively. */
 void Gauss_Quadrature::open_Interval( const double x_zero, const char openside )
 {
    cout<<"  Gauss_tag0"<<x_zero<<endl;
    unit(T);
    cout<<"  Gauss_tag1"<<x_zero<<endl;
    switch(openside)
    {
        case 'L':
            std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) \
                 << "    Transforms to (infty," << x_zero << "].\n\n";
            transGaussq_unit2inf_x0( x_zero );
            break;
        case 'R':
            cout<<"  Gauss_tag2"<<x_zero<<endl;
            //std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(2) \
                 << "    Transforms to [" << x_zero << ", infty)\n\n";
            cout<< "    Transforms to [" << x_zero << ", infty)\n\n";
            cout<<"  Gauss_tag1"<<x_zero<<endl;
            transGaussq_unit2x0_inf( x_zero );
        break;
        default:
        std::cout << "\n  ! Wrong assignment when indicating the side of the half-open interval in which the infinity " \
                    "is set!\n";
        std::cout << "\n    Available choices:\n\t\'L\': Open on the left-hand side of the interval i.e., starts with " \
                    "negative infinity.\n\t\'R\': Open on the right-hand side of the interval i.e., ended with "       \
                    "positive infinity.";
        break;
    }

   return;
 }

/*Contructor.
   --Generates T-type Gauss quadrature in the interval (-infty,infty). */
 void Gauss_Quadrature::infty_Interval( const char* inf_range )
 {
   unit(T);
   if(  inf_range == "inf" || inf_range =="infty" )
   {
     std::cout << "    Transforms to (-infty, infty)\n\n";
     transGaussq_unit2inf();
   }
   else if( inf_range == "0-inf" || inf_range == "0-infty" )
   {
     std::cout << "    Transforms to [0.0,infty)\n\n";
     transGaussq_unit20_inf();
   }
   else
     std::cout << "\n  ! Wrong assignment for the open interval of infinity!\n\n Try \"inf\" or \"0-inf\".\n\n";

   return;
 }

/* Destructor. */
 Gauss_Quadrature::~Gauss_Quadrature()
 {
   free();
 }

/* Free data members' workspace. */
 void Gauss_Quadrature::free()
 {
   if( N > 0 )
     delete [] nodes, weights;

   return;
 }

/* * * Implementations * * * * * * * * * * * * * * * * * * * * * * * * */

/*PURPOSE
   --Generate the Legendre polynomials and the corresponding derivatives.
     The weight function w(x) = 1.                                     */
 void Gauss_Quadrature::Legendre_polynomial( const int N, const double x, double &y, double &dy )
 {
   double p0, p, dp, rtmp, xtmp;
   int n;

   p0 = 1.0;    p = x;
   xtmp = 1.0/ ( x*x -1.0 );
   /* Recurrence relation:  for n = 1, 2, ..., N,
   *    nL_n(x)  = (2n-1)xL_{n-1}(x) - (n-1)L_{n-2}(x)
   *    (x^2-1)L'_n(x) = nxL_n(x) - nL_{n-1}(x)
   */
   for( n = 1; n != N; ++n )
   {
     rtmp = p0;    p0 = p;
      p = ( ( 2*n +1 ) *x *p0 - n *rtmp ) /(n+1);
     dp = (n+1) * xtmp * ( x*p - p0 );
   }
   y = p; dy = dp;

   return;
 }

/*PURPOSE
   --Generate the Gauss-Legendre Quadrature on the interval [-1,1].
       N
        --The number of Gauss-Legendre grids.
       Gq
        --Keep the Gauss-Legendre roots and corresponding weights (Gq[0][:], Gq[1][:]).   */
 void Gauss_Quadrature::Gaussq_Legendre()
 {
   double x, y, dy, ratio, N_;
   int i;

   N_ = ( N +1 ) /2;
   for( i = 0; i != N_; ++i )
   {
/*    Using initial guess: x_zero = cos( PI* (i-0.25)/(N+0.5) ), i = 1, 2, ..., N. OR: x_zero = cos( PI* (k+0.75)/(N+0.5) ) */
      x = cos( PI *( i+1 -0.25 ) / ( N +0.5 ) );
      
/*    Search the root by Newton-Raphson method. */
      do
      {
        Legendre_polynomial( N, x, y, dy );
        ratio = y/dy;
        x = x - ratio;
      } while( fabs(ratio) > 1e-16 );
      nodes[i] = -x;
      weights[i] = 2.0/( ( 1.0-x*x) * dy*dy );
      
/*    Applying the symmetry of Gauss-Legendre nodes to reduce computation. */
      nodes[N-1-i] = -nodes[i];
      weights[N-1-i] =  weights[i];
    }

   return;
 }

/*PURPOSE
   --Generate the Gauss-Chebyshev Quadrature based on Chebyshev points of the first type (kind).
     The weght function w(x) = 1/ sqrt(1-x^2).
     The zeros of T_N(x) in [-1,1].
     ix = cos( (i+1/2)PI /N ), iwx = PI /N with, i = 0, 1, ..., N-1.                          */
 void Gauss_Quadrature::Gaussq_ChebyshevT1()
 {
   int i, M;

   M = (N+1)/2;
   for( i = 0; i != M; ++i )
   {
      nodes[i] = cos( PI*  ( i +0.5 ) /N );
      weights[i] = sqrt( 1.0- nodes[i]*nodes[i] ) *PI /N;
/*    Apply the symmetry property of the nodes. */
      nodes[N-1-i] = -nodes[i];
      weights[N-1-i] =  weights[i];
   }

   return;
 }

/*PURPOSE
   --Get the Gauss-Chebyshev Quadrature based on Chebyshev points of the second  type (kind), i.e., the nodes are the
     local extrema of T_N(x) in [-1,1].
       ix = cos( i*PI /(N+1) )
       iwx = PI/(N+1) *sin^2 ( i*PI /(N+1) ), with i = 1, ..., N.                                                  */
 void Gauss_Quadrature::Gaussq_ChebyshevT2()
 {
   int i, M, N_;

   M = (N +1)/2;    N_ = N +1;
   for( i = 0; i != M; ++i )
   {
      nodes[i] = cos( ( i +1 ) *PI /N_ );
      weights[i] = 1.0/ sqrt( 1.0- nodes[i]*nodes[i] ) *PI /N_ \
                   *pow( sin( PI* ( i +1.0 ) /N_ ), 2 );
      nodes[N-1-i] = -nodes[i];
      weights[N-1-i] =  weights[i];
   }

   return;
 }

/*PURPOSE
   --Linearly transforming Gauss Quadrature from [-1,1] to [x_zero,x_one].
        x' = (b-a)/2 *x + (a+b)/2
       wx' = (b-a)/2 *wx                   */
 void Gauss_Quadrature::transGaussq_linear( const double x_zero, const double x_one )
 {
   int i;
   double c1, c2;

   c1 = 0.5* (x_one+x_zero);   c2 = 0.5* (x_one-x_zero);
   for( i = 0; i != N; ++i )
   {
      nodes[i] = c2* nodes[i] +c1;
      weights[i] = c2* weights[i];
   }

   return;
 }

/*PURPOSE
   --Transforming the Gauss-Legendre Quadrature from [-1,1] to [0,inf). The modification proposed for the evaluation of
     the Casimir-Polder integral is applied.
       wx' = 1/(1-x)^2 *wx
        x' = 0.5* (1+x)/(1-x)                  */
 void Gauss_Quadrature::transGaussq_unit20_inf()
 {
   double c0;
   int i;

   c0 = 0.5;
   for( i = 0; i != N; ++i )
   {
      weights[i] = 2.0* c0* weights[i] / ( ( 1.0- nodes[i] ) * ( 1.0- nodes[i] ) );
      nodes[i] = c0* ( 1.0+ nodes[i] ) / ( 1.0- nodes[i] );
   }

   return;
 }

/*PURPOSE                                                               
   --Transforming Gauss Quadrature from [-1,1] to (-inf,inf).
       wx' = 5* (1+x^2)/(1-x^2)^2 *wx
        x' = 5* x/(1-x^2)                                  */
 void Gauss_Quadrature::transGaussq_unit2inf()
 {
   int i;
   double xsq;

   for( i = 0; i != N; ++i )
   {
      xsq = nodes[i]*nodes[i];
      weights[i] = weights[i] *5.0 *( 1.0+ xsq ) /( (1.0- xsq)*(1.0- xsq) );
      nodes[i] = 5.0* nodes[i] / ( 1.0- xsq );
   }

   return;
 }

/*PURPOSE
   --Transforming Gauss Quadrature from [-1,1] to [x_zero,inf).
       wx' = 30/(1-x)^2 *wx
        x' = 15* (x+1)/(1-x) +x_zero                         */
 void Gauss_Quadrature::transGaussq_unit2x0_inf( const double x_zero )
 {
   int i;

   for( i = 0; i != N; ++i )
   {
/*    weights[i] = weights[i] *30.0 /( ( 1.0- nodes[i] )*( 1.0- nodes[i] ) );
      nodes[i] = 15.0* ( nodes[i] +1.0 ) / ( 1.0- nodes[i] ) + x_zero;         */
      weights[i] = weights[i] /( ( 1.0- nodes[i] )*( 1.0- nodes[i] ) );
      nodes[i] = 0.5* ( nodes[i] +1.0 ) / ( 1.0- nodes[i] ) + x_zero;
   }

   return;
 }

/*PURPOSE
   --Transforming Gauss Quadrature from [-1,1] to (-inf,x_zero].
       wx' = 30 /(1+x)^2 *wx 
        x' = 15* (x-1)/(1+x) +x_zero                          */
 void Gauss_Quadrature::transGaussq_unit2inf_x0( const double x_zero )
 {
   int i;

   for( i = 0; i != N; ++i )
   {
/*    weights[i] = weights[i] *30.0 /(  ( 1.0+ nodes[i] )*( 1.0+ nodes[i] ) );
      nodes[i] = 15.0* ( nodes[i] -1.0 ) / ( 1.0+ nodes[i] ) + x_zero;            */
      weights[i] = weights[i] /(  ( 1.0+ nodes[i] )*( 1.0+ nodes[i] ) );
      nodes[i] = 0.5* ( nodes[i] -1.0 ) / ( 1.0+ nodes[i] ) + x_zero;
   }

   return;
 }

/*PURPOSE
   --Transforming the Gauss Quadrature from (-inf,inf) to [-1,1].
       wx' = wx* 10* ( 1+ 1/ sqrt( 25+ 4x^2 ) ) /( 5+ sqrt( 25+ 4x^2 ) )^2
        x' = 2x/ ( 5+ sqrt( 25+ 4x^2 ) )                                  */
 void Gauss_Quadrature::transGaussq_inf2unit()
 {
   int i;
   double rtmp, xsq;

   for( i = 0; i != N; ++i )
   {
      xsq = 4.0* nodes[i]*nodes[i];
      rtmp = sqrt( 25.0+ xsq );
      weights[i] = weights[i]* 10.0* ( 1.0+ 1.0/rtmp ) /( ( 5.0+ rtmp ) *( 5.0+ rtmp ) );
      nodes[i] = 2.0*nodes[i] /( 5.0+ sqrt( 25.0+ xsq ) );
   }

   return;
 }

/*PURPOSE
   --Transforming the Gauss Quadrature from [x_zero,inf) to [-1,1].
       wx' = wx* 30/(x+15-x_zero)^2
        x' = (x-15-x_zero)/(x+15-x_zero)                             */
 void Gauss_Quadrature::transGaussq_x0_inf2unit( const double x_zero )
 {
   int i;
   double dtmp;

   for( i = 0; i != N; ++i )
   {
      dtmp = 1.0/ ( nodes[i] - 15.0 -x_zero);
      weights[i] = weights[i] *30.0 *dtmp *dtmp;
      nodes[i] = ( nodes[i] - 15.0 -x_zero ) *dtmp;
   }

   return;
 }

/*PURPOSE
   --Transforming the Gauss Quadrature from (-inf,x_zero] to [-1,1].
     wx' = wx* 30/(15-x+x_zero)^2
      x' = (15+x-x_zero)/(15-x+x_zero)                                */
 void Gauss_Quadrature::transGaussq_inf_x02unit( const double x_zero )
 {
   int i;
   double dtmp;

   for( i = 0; i != N; ++i )
   {
      dtmp = 1.0/ ( 15.0- nodes[i] +x_zero );
      weights[i] = weights[i] *30.0 *dtmp *dtmp;
      nodes[i] = ( 15.0+ nodes[i] -x_zero ) *dtmp;
   }

   return;
 }
