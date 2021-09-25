#ifndef INTERPOLATION_H
#define INTERPOLATION_H

/**
 @file   interpolation.h
 @author James Gregson james.gregson@gmail.com
 @date   09/03/2012
 
 Simple templated implementation of several classes of interpolation routines in 
 one, two and three dimensions.  Each class provides several options regarding 
 input type and may optionally return a vector of weights associated with each
 point.
 
 Free for all use with proper attribution, e.g. by keeping this header intact
 for source distributions or by brief mention in binary distribution documentation.
 Please contact me if you find this useful, would like features added and/or 
 find bugs or want to contribute fixes.
 
 Interpolation types currently implemented are:
 
	- Linear (*_lerp_* functions) 1/2/3D functions with higher dimensions
      obtained via tensor product
 
    - Cosine (*_cosine_* functions) 1/2/3D functions with higher dimensions
      obtained via tensor product. Produces smooth interpolants but cannot
      reproduce linear gradients.
 
    - Cubic (*_cubic_* functions) 1/2/3D functions with higher dimensions
      obtains via tensor product. More limited API in 3D due to large number 
      of input values required (64 for 3D). Uses Catmull-Rom weights to
      give C1 continuity of adjacent segments.
 
 All methods use a common notation for input values, v000 is at 
 [x,y,z] = [0,0,0], while v123 is at [x,y,z] = [1, 2, 3]. Parameters packed 
 into arrays are ordered by scanline, then slice, e.g. for a 3D voxel the packed
 array of corners would be:
 
 v = [ v000, v100, v010, v110, v001, v101, v011, v111 ]
 
 All interpolation methods are tested for consistency between the various
 API functions, and that they reproduce their input values properly for
 appropriate parametric values.  The test suite may be run by calling 
 the function interpolation_test()
*/

// needed for cos() and M_PI
#include<cmath>

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 1D Linear Interpolation ////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes weights for 1D linear interpolation
 @param[in]  t  fractional distance along interval
 @param[out] w0 left value
 @param[out] w1 right value
*/
template< typename S >
inline void interpolation_lerp_1D_weights( const S &t, S &w0, S &w1 ){
	w0 = S(1.0)-t;
	w1 = t;
}

/**
 @brief performs 1D linear interpolation using weights computed by interpolation_lerp_1D_weights()
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_1D( const V &v0, const V &v1, const S &t ){
	S w0, w1;
	interpolation_lerp_1D_weights( t, w0, w1 );
	return v0*w0 + v1*w1;
}

/**
 @brief performs 1D linear interpolation using weights computed by interpolation_lerp_1D_weights(),
 also returning the weights that were used
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @param[out] w0 left value weight
 @param[out] w1 right value weight
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_1D( const V &v0, const V &v1, const S &t, S &w0, S &w1 ){
	interpolation_lerp_1D_weights( t, w0, w1 );
	return v0*w0 + v1*w1;
}

/**
 @brief performs 1D linear interpolation using weights computed by interpolation_lerp_1D_weights()
 taking input values as an array
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_1D( const V *v, const S &t ){
	S w0, w1;
	interpolation_lerp_1D_weights( t, w0, w1 );
	return v[0]*w0 + v[1]*w1;
}

/**
 @brief performs 1D linear interpolation using weights computed by interpolation_lerp_1D_weights()
 returning weights, taking inputs and outputs as arrays
 @param[in]  v  array of input values as v = [ left_val, right_val ]
 @param[in]  t  fractional distance along interval
 @param[out] w  array of output weights, w = [ w(left_val), w(right_val) ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_1D( const V *v, const S &t, S *w ){
	interpolation_lerp_1D_weights( t, w[0], w[1] );
	return v[0]*w[0] + v[1]*w[1];
}

/**
 @brief test function for the 1D linear interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input, and that they reproduce the input values for 
 appropriate parametric coordinates
 @return true if test passes false otherwise
*/
static bool interpolation_test_lerp_1D(){
	const double t=0.325f, v[] = { 1.0, 2.0 };
	double o0, o1, o2, o3, w[2];
	bool result = true;
	
	// consistency check between different API functions
	o0 = interpolation_lerp_1D( v[0], v[1], t );
	o1 = interpolation_lerp_1D( v[0], v[1], t, w[0], w[1] );
	o2 = interpolation_lerp_1D( v, t );
	o3 = interpolation_lerp_1D( v, t, w );
	result = fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check that left value is reproduced
	o0 = interpolation_lerp_1D( v[0], v[1], 0.0 );
	result &= fabs(o0-v[0]) < 1e-9;
	
	// check that right value is reproduced
	o0 = interpolation_lerp_1D( v[0], v[1], 1.0 );
	result &= fabs(o0-v[1]) < 1e-9;
	
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 1D Cosine Interpolation (Kind of useless) //////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 @brief computes weights for 1D cosine interpolation
 @param[in]  t  fractional distance along interval
 @param[out] w0 left value
 @param[out] w1 right value
*/
template< typename S >
inline void interpolation_cosine_1D_weights( const S &t, S &w0, S &w1 ){
	w0 = S(cos(double(t)*M_PI)/2.0+0.5);
	w1 = S(1.0)-w0;
}

/**
 @brief performs 1D cosine interpolation using weights computed by interpolation_cosine_1D_weights()
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_1D( const V &v0, const V &v1, const S &t ){
	S w0, w1;
	interpolation_cosine_1D_weights( t, w0, w1 );
	return v0*w0 + v1*w1;
}

/**
 @brief performs 1D cosine interpolation using weights computed by interpolation_cosine_1D_weights(),
 also returning the weights that were used
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @param[out] w0 left value weight
 @param[out] w1 right value weight
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_1D( const V &v0, const V &v1, const S &t, S &w0, S &w1 ){
	interpolation_cosine_1D_weights( t, w0, w1 );
	return v0*w0 + v1*w1;
}

/**
 @brief performs 1D cosine interpolation using weights computed by interpolation_cosine_1D_weights()
 taking input values as an array
 @param[in]  v  input values v=[v0, v1]
 @param[in]  t  fractional distance along interval
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_1D( const V *v, const S &t ){
	S w0, w1;
	interpolation_cosine_1D_weights( t, w0, w1 );
	return v[0]*w0 + v[1]*w1;
}

/**
 @brief performs 1D cosine interpolation using weights computed by interpolation_cosine_1D_weights()
 returning weights, taking inputs and outputs as arrays
 @param[in]  v  array of input values as v = [ left_val, right_val ]
 @param[in]  t  fractional distance along interval
 @param[out] w  array of output weights, w = [ w(left_val), w(right_val) ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_1D( const V *v, const S &t, S *w ){
	interpolation_cosine_1D_weights( t, w[0], w[1] );
	return v[0]*w[0] + v[1]*w[1];
}

/**
 @brief test function for the 1D cosine interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
*/
static bool interpolation_test_cosine_1D(){
	const double t=0.325f, v[] = { 1.0, 2.0 };
	double o0, o1, o2, o3, w[2];
	bool result = true;
	
	// test that API functions are consistent
	o0 = interpolation_cosine_1D( v[0], v[1], t );
	o1 = interpolation_cosine_1D( v[0], v[1], t, w[0], w[1] );
	o2 = interpolation_cosine_1D( v, t );
	o3 = interpolation_cosine_1D( v, t, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check that left value is reproduced
	o0 = interpolation_cosine_1D( v[0], v[1], 0.0 );
	result &= fabs(o0-v[0]) < 1e-9;
	
	// check that right value is reproduced
	o0 = interpolation_cosine_1D( v[0], v[1], 1.0 );
	result &= fabs(o0-v[1]) < 1e-9;
	
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 1D Cubic Interpolation (Catmull-Rom Weights) ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////


/**
 @brief computes weights for 1D cubic Catmull-Rom spline interpolation
 @param[in]  t  fractional distance along interval
 @param[out] 
*/
template< typename S >
inline void interpolation_cubic_1D_weights( const S &t, S &w0, S &w1, S &w2, S &w3 ){
	S t2=t*t, t3=t2*t;
	// Catmull-Rom spline weights
	//    cubic       quadratic   linear      constant
	w0 = -S(0.5)*t3 + S(1.0)*t2 - S(0.5)*t;
	w1 =  S(1.5)*t3 - S(2.5)*t2             + S(1.0);
	w2 = -S(1.5)*t3 + S(2.0)*t2 + S(0.5)*t;
	w3 =  S(0.5)*t3 - S(0.5)*t2;
}

/**
 @brief performs 1D cosine interpolation using weights computed by interpolation_cosine_1D_weights()
 @param[in]  v0 leftmost value
 @param[in]  v1 immediate left value
 @param[in]  v2 immediate right value
 @param[in]  v3 rightmost value
 @param[in]  t  fractional distance along interval
 @return interpolated value
*/
template< typename V, typename S >
inline V interpolation_cubic_1D( const V &v0, const V &v1, const V &v2, const V &v3, const S &t ){
	S w0, w1, w2, w3;
	interpolation_cubic_1D_weights( t, w0, w1, w2, w3 );
	return v0*w0 + v1*w1 + v2*w2 + v3*w3;
}

/**
 @brief performs 1D cubic interpolation using weights computed by interpolation_cubic_1D_weights(),
 also returning the weights that were used
 @param[in]  v0 left value
 @param[in]  v1 right value
 @param[in]  t  fractional distance along interval
 @param[out] w0 left value weight
 @param[out] w1 right value weight
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cubic_1D( const V &v0, const V &v1, const V &v2, const V &v3, const S &t, S &w0, S &w1, S &w2, S &w3 ){
	interpolation_cubic_1D_weights( t, w0, w1, w2, w3 );
	return v0*w0 + v1*w1 + v2*w2 + v3*w3;
}

/**
 @brief performs 1D cubic interpolation using weights computed by interpolation_cubic_1D_weights()
 taking input values as an array
 @param[in]  v  input values v = [v0, v1, v2, v3]
 @param[in]  t  fractional distance along interval
 @return interpolated value
*/
template< typename V, typename S >
inline V interpolation_cubic_1D( const V *v, const S &t ){
	S w0, w1, w2, w3;
	interpolation_cubic_1D_weights( t, w0, w1, w2, w3 );
	return v[0]*w0 + v[1]*w1 + v[2]*w2 + v[3]*w3;
}

/**
 @brief performs 1D cubic interpolation using weights computed by interpolation_cubic_1D_weights()
 returning weights, taking inputs and outputs as arrays
 @param[in]  v  array of input values as v = [ left_val, right_val ]
 @param[in]  t  fractional distance along interval
 @param[out] w  array of output weights, w = [ w(left_val), w(right_val) ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cubic_1D( const V *v, const S &t, S *w ){
	interpolation_cubic_1D_weights( t, w[0], w[1], w[2], w[3] );
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] + v[3]*w[3];
}

/**
 @brief test function for the 1D cubic interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_cubic_1D(){
	const double t=0.325f, v[] = { 2.0, 1.0, 2.0, -1.0 };
	double o0, o1, o2, o3, w[4];
	bool result = true;
	
	// check that API functions are consistent
	o0 = interpolation_cubic_1D( v[0], v[1], v[2], v[3], t );
	o1 = interpolation_cubic_1D( v[0], v[1], v[2], v[3], t, w[0], w[1], w[2], w[3] );
	o2 = interpolation_cubic_1D( v, t );
	o3 = interpolation_cubic_1D( v, t, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check that left value is reproduced
	o0 = interpolation_cubic_1D( v[0], v[1], v[2], v[3], 0.0 );
	result &= fabs(o0-v[1]) < 1e-9;
	
	// check that right value is reproduced
	o0 = interpolation_cubic_1D( v[0], v[1], v[2], v[3], 1.0 );
	result &= fabs(o0-v[2]) < 1e-9;
	
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 2D Linear Interpolation (Tensor product formulation) ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for 2D bilinear interpolation
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1]
*/
template< typename S >
inline void interpolation_lerp_2D_weights( const S &t0, const S &t1, S &w00, S &w10, S &w01, S &w11 ){
	S wlo, whi;
	interpolation_lerp_1D_weights( t0, w00, w10 );
	interpolation_lerp_1D_weights( t1, wlo, whi );
	w01=w00; w11=w10;
	w00*=wlo; w10*=wlo;
	w01*=whi; w11*=whi;
}

/**
 @brief perform 2D linear interpolation using weights computed by interpolation_lerp_2D_weights()
 @param[in] v00 value for [x,y] = [0,0]
 @param[in] v10 value for [x,y] = [1,0]
 @param[in] v01 value for [x,y] = [0,1]
 @param[in] v11 value for [x,y] = [1,1]
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_2D( const V &v00, const V &v10, const V &v01, const V &v11, const S &t0, const S &t1 ){
	S w00, w10, w01, w11;
	interpolation_lerp_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v00*w00 + v10*w10 + v01*w01 + v11*w11;
}

/**
 @brief perform 2D linear interpolation using weights computed by interpolation_lerp_2D_weights()
 and returns the weights that are used
 @param[in]  v00 value for [x,y] = [0,0]
 @param[in]  v10 value for [x,y] = [1,0]
 @param[in]  v01 value for [x,y] = [0,1]
 @param[in]  v11 value for [x,y] = [1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1] 
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_2D( const V &v00, const V &v10, const V &v01, const V &v11, const S &t0, const S &t1, S &w00, S &w10, S& w01, S &w11 ){
	interpolation_lerp_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v00*w00 + v10*w10 + v01*w01 + v11*w11;
}

/**
 @brief perform 2D linear interpolation using weights computed by interpolation_lerp_2D_weights(),
 taking the input parameters as a vector
 @param[in] v   array of input values v = [ v00, v10, v01, v11 ]
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_2D( const V *v, const S &t0, const S &t1 ){
	S w00, w10, w01, w11;
	interpolation_lerp_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v[0]*w00 + v[1]*w10 + v[2]*w01 + v[3]*w11;
}

/**
 @brief perform 2D linear interpolation using weights computed by interpolation_lerp_2D_weights() and
 returns weights, taking the input and output parameters as arrays
 @param[in]  v   array of input values v = [ v00, v10, v01, v11 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w   array of weight values w = [ w00, w10, w01, w11 ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_2D( const V *v, const S &t0, const S &t1, S *w ){
	interpolation_lerp_2D_weights( t0, t1, w[0], w[1], w[2], w[3] );
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] + v[3]*w[3];
}

/**
 @brief test function for the 2D linear interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
*/
static bool interpolation_test_lerp_2D(){
	const double t0=0.325f, t1=0.675, v[] = { 1.0, 2.0, 2.0, 5.0 };
	double o0, o1, o2, o3, w[4];
	bool result = true;
	
	// check that API functions are consistent
	o0 = interpolation_lerp_2D( v[0], v[1], v[2], v[3], t0, t1 );
	o1 = interpolation_lerp_2D( v[0], v[1], v[2], v[3], t0, t1, w[0], w[1], w[2], w[3] );
	o2 = interpolation_lerp_2D( v, t0, t1 );
	o3 = interpolation_lerp_2D( v, t0, t1, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check that input values are reproduced
	result &= fabs( interpolation_lerp_2D( v[0], v[1], v[2], v[3], 0.0, 0.0 )-v[0] ) < 1e-9; 
	result &= fabs( interpolation_lerp_2D( v[0], v[1], v[2], v[3], 1.0, 0.0 )-v[1] ) < 1e-9;
	result &= fabs( interpolation_lerp_2D( v[0], v[1], v[2], v[3], 0.0, 1.0 )-v[2] ) < 1e-9;
	result &= fabs( interpolation_lerp_2D( v[0], v[1], v[2], v[3], 1.0, 1.0 )-v[3] ) < 1e-9; 

	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 2D Cosine Interpolation (Tensor product formulation) ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for 2D cosine interpolation
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1]
 */
template< typename S >
inline void interpolation_cosine_2D_weights( const S &t0, const S &t1, S &w00, S &w10, S &w01, S &w11 ){
	S wlo, whi;
	interpolation_cosine_1D_weights( t0, w00, w10 );
	interpolation_cosine_1D_weights( t1, wlo, whi );
	w01=w00; w11=w10;
	w00*=wlo; w10*=wlo;
	w01*=whi; w11*=whi;
}

/**
 @brief perform 2D cosine interpolation using weights computed by interpolation_cosine_2D_weights()
 @param[in] v00 value for [x,y] = [0,0]
 @param[in] v10 value for [x,y] = [1,0]
 @param[in] v01 value for [x,y] = [0,1]
 @param[in] v11 value for [x,y] = [1,1]
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_2D( const V &v00, const V &v10, const V &v01, const V &v11, const S &t0, const S &t1 ){
	S w00, w10, w01, w11;
	interpolation_cosine_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v00*w00 + v10*w10 + v01*w01 + v11*w11;
}

/**
 @brief perform 2D cosine interpolation using weights computed by interpolation_cosine_2D_weights()
 and returns the weights that are used
 @param[in]  v00 value for [x,y] = [0,0]
 @param[in]  v10 value for [x,y] = [1,0]
 @param[in]  v01 value for [x,y] = [0,1]
 @param[in]  v11 value for [x,y] = [1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1] 
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_2D( const V &v00, const V &v10, const V &v01, const V &v11, const S &t0, const S &t1, S &w00, S &w10, S& w01, S &w11 ){
	interpolation_cosine_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v00*w00 + v10*w10 + v01*w01 + v11*w11;
}

/**
 @brief perform 2D cosine interpolation using weights computed by interpolation_cosine_2D_weights(),
 taking the input parameters as a vector
 @param[in] v   array of input values v = [ v00, v10, v01, v11 ]
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_2D( const V *v, const S &t0, const S &t1 ){
	S w00, w10, w01, w11;
	interpolation_cosine_2D_weights( t0, t1, w00, w10, w01, w11 );
	return v[0]*w00 + v[1]*w10 + v[2]*w01 + v[3]*w11;
}

/**
 @brief perform 2D cosine interpolation using weights computed by interpolation_cosine_2D_weights() and
 returns weights, taking the input and output parameters as arrays
 @param[in]  v   array of input values v = [ v00, v10, v01, v11 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w   array of weight values w = [ w00, w10, w01, w11 ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_2D( const V *v, const S &t0, const S &t1, S *w ){
	interpolation_cosine_2D_weights( t0, t1, w[0], w[1], w[2], w[3] );
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] + v[3]*w[3];
}

/**
 @brief test function for the 2D cosine interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_cosine_2D(){
	const double t0=0.325f, t1=0.675, v[] = { 1.0, 2.0, 2.0, 5.0 };
	double o0, o1, o2, o3, w[4];
	bool result = true;
	
	// check API consistency
	o0 = interpolation_cosine_2D( v[0], v[1], v[2], v[3], t0, t1 );
	o1 = interpolation_cosine_2D( v[0], v[1], v[2], v[3], t0, t1, w[0], w[1], w[2], w[3] );
	o2 = interpolation_cosine_2D( v, t0, t1 );
	o3 = interpolation_cosine_2D( v, t0, t1, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check that input values are reproduced
	result &= fabs( interpolation_cosine_2D( v[0], v[1], v[2], v[3], 0.0, 0.0 )-v[0] ) < 1e-9; 
	result &= fabs( interpolation_cosine_2D( v[0], v[1], v[2], v[3], 1.0, 0.0 )-v[1] ) < 1e-9;
	result &= fabs( interpolation_cosine_2D( v[0], v[1], v[2], v[3], 0.0, 1.0 )-v[2] ) < 1e-9;
	result &= fabs( interpolation_cosine_2D( v[0], v[1], v[2], v[3], 1.0, 1.0 )-v[3] ) < 1e-9; 
	
	return result;
	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 2D Cubic Interpolation (Tensor product formulation) ////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for 2D cubic interpolation
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w20 weight for [x,y] = [2,0]
 @param[out] w30 weight for [x,y] = [3,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1]
 @param[out] w21 weight for [x,y] = [2,1]
 @param[out] w31 weight for [x,y] = [3,1]
 @param[out] w02 weight for [x,y] = [0,2]
 @param[out] w12 weight for [x,y] = [1,2]
 @param[out] w22 weight for [x,y] = [2,2]
 @param[out] w32 weight for [x,y] = [3,2]
 @param[out] w03 weight for [x,y] = [0,3]
 @param[out] w13 weight for [x,y] = [1,3]
 @param[out] w23 weight for [x,y] = [2,3]
 @param[out] w33 weight for [x,y] = [3,3] 
*/
template< typename S >
inline void interpolation_cubic_2D_weights( const S &t0, const S &t1, S &w00, S &w10, S &w20, S &w30, S &w01, S &w11, S &w21, S &w31, S &w02, S &w12, S &w22, S &w32, S &w03, S &w13, S &w23, S &w33 ){
	S tw0, tw1, tw2, tw3;
	interpolation_cubic_1D_weights( t0, w00, w10, w20, w30 );
	interpolation_cubic_1D_weights( t1, tw0, tw1, tw2, tw3 );
	w01=w00; w11=w10; w21=w20; w31=w30;
	w02=w00; w12=w10; w22=w20; w32=w30;
	w03=w00; w13=w10; w23=w20; w33=w30;	
	w00*=tw0; w10*=tw0; w20*=tw0; w30*=tw0;
	w01*=tw1; w11*=tw1; w21*=tw1; w31*=tw1;
	w02*=tw2; w12*=tw2; w22*=tw2; w32*=tw2;
	w03*=tw3; w13*=tw3; w23*=tw3; w33*=tw3;
}

/**
 @brief perform 2D cubic interpolation using weights computed by interpolation_cubic_2D_weights()
 @param[in] v00 value for [x,y] = [0,0]
 @param[in] v10 value for [x,y] = [1,0]
 @param[in] v20 value for [x,y] = [2,0]
 @param[in] v30 value for [x,y] = [3,0]
 @param[in] v01 value for [x,y] = [0,1]
 @param[in] v11 value for [x,y] = [1,1]
 @param[in] v21 value for [x,y] = [2,1]
 @param[in] v31 value for [x,y] = [3,1]
 @param[in] v02 value for [x,y] = [0,2]
 @param[in] v12 value for [x,y] = [1,2]
 @param[in] v22 value for [x,y] = [2,2]
 @param[in] v32 value for [x,y] = [3,2]
 @param[in] v02 value for [x,y] = [0,2]
 @param[in] v12 value for [x,y] = [1,2]
 @param[in] v22 value for [x,y] = [2,2]
 @param[in] v32 value for [x,y] = [3,2] 
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
*/
template< typename V, typename S >
inline V interpolation_cubic_2D( const V &v00, const V &v10, const V &v20, const V &v30, const V &v01, const V &v11, const V &v21, const V &v31, const V &v02, const V &v12, const V &v22, const V &v32, const V &v03, const V &v13, const V &v23, const V &v33, const S &t0, const S &t1 ){
	S w00, w10, w20, w30, w01, w11, w21, w31, w02, w12, w22, w32, w03, w13, w23, w33;
	interpolation_cubic_2D_weights( t0, t1, w00, w10, w20, w30, w01, w11, w21, w31, w02, w12, w22, w32, w03, w13, w23, w33 );
	return v00*w00 + v10*w10 + v20*w20 + v30*w30
	     + v01*w01 + v11*w11 + v21*w21 + v31*w31
	     + v02*w02 + v12*w12 + v22*w22 + v32*w32
	     + v03*w03 + v13*w13 + v23*w23 + v33*w33;
}

/**
 @brief perform 2D cubic interpolation using weights computed by interpolation_cubic_2D_weights()
 @param[in] v00 value for [x,y] = [0,0]
 @param[in] v10 value for [x,y] = [1,0]
 @param[in] v20 value for [x,y] = [2,0]
 @param[in] v30 value for [x,y] = [3,0]
 @param[in] v01 value for [x,y] = [0,1]
 @param[in] v11 value for [x,y] = [1,1]
 @param[in] v21 value for [x,y] = [2,1]
 @param[in] v31 value for [x,y] = [3,1]
 @param[in] v02 value for [x,y] = [0,2]
 @param[in] v12 value for [x,y] = [1,2]
 @param[in] v22 value for [x,y] = [2,2]
 @param[in] v32 value for [x,y] = [3,2]
 @param[in] v02 value for [x,y] = [0,2]
 @param[in] v12 value for [x,y] = [1,2]
 @param[in] v22 value for [x,y] = [2,2]
 @param[in] v32 value for [x,y] = [3,2] 
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @param[out] w00 weight for [x,y] = [0,0]
 @param[out] w10 weight for [x,y] = [1,0]
 @param[out] w20 weight for [x,y] = [2,0]
 @param[out] w30 weight for [x,y] = [3,0]
 @param[out] w01 weight for [x,y] = [0,1]
 @param[out] w11 weight for [x,y] = [1,1]
 @param[out] w21 weight for [x,y] = [2,1]
 @param[out] w31 weight for [x,y] = [3,1]
 @param[out] w02 weight for [x,y] = [0,2]
 @param[out] w12 weight for [x,y] = [1,2]
 @param[out] w22 weight for [x,y] = [2,2]
 @param[out] w32 weight for [x,y] = [3,2]
 @param[out] w02 weight for [x,y] = [0,2]
 @param[out] w12 weight for [x,y] = [1,2]
 @param[out] w22 weight for [x,y] = [2,2]
 @param[out] w32 weight for [x,y] = [3,2]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cubic_2D( const V &v00, const V &v10, const V &v20, const V &v30, const V &v01, const V &v11, const V &v21, const V &v31, const V &v02, const V &v12, const V &v22, const V &v32, const V &v03, const V &v13, const V &v23, const V &v33, const S &t0, const S &t1, S &w00, S &w10, S &w20, S &w30, S &w01, S &w11, S &w21, S &w31, S &w02, S &w12, S &w22, S &w32, S &w03, S &w13, S &w23, S &w33 ){
	interpolation_cubic_2D_weights( t0, t1, w00, w10, w20, w30, w01, w11, w21, w31, w02, w12, w22, w32, w03, w13, w23, w33 );
	return v00*w00 + v10*w10 + v20*w20 + v30*w30
	+ v01*w01 + v11*w11 + v21*w21 + v31*w31
	+ v02*w02 + v12*w12 + v22*w22 + v32*w32
	+ v03*w03 + v13*w13 + v23*w23 + v33*w33;
}

/**
 @brief perform 2D cubic interpolation using weights computed by interpolation_cubic_2D_weights(),
 taking the input parameters as a vector
 @param[in] v   array of input values v = [ v00, v10, v20, v30, v01, v11, v21, v31, v02, v12, v22, v32, v03, v13, v23, v33 ]
 @param[in] t0  fractional distance along first axis
 @param[in] t1  fractional distance along second axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cubic_2D( const V *v, const S &t0, const S &t1 ){
	S sum=S(0.0), w[16];
	interpolation_cubic_2D_weights( t0, t1, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9], w[10], w[11], w[12], w[13], w[14], w[15] );
	for( int i=0; i<16; i+=4 ){
		sum += v[i+0]*w[i+0];
		sum += v[i+1]*w[i+1];
		sum += v[i+2]*w[i+2];
		sum += v[i+3]*w[i+3];
	}
	return sum;
}

/**
 @brief perform 2D cubic interpolation using weights computed by interpolation_cubic_2D_weights(),
 taking the input parameters as a vector
 @param[in]  v   array of input values v = [ v00, v10, v20, v30, v01, v11, v21, v31, v02, v12, v22, v32, v03, v13, v23, v33 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[out] w   array of output weights w = [ w00, w10, w20, w30, w01, w11, w21, w31, w02, w12, w22, w32, w03, w13, w23, w33 ]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cubic_2D( const V *v, const S &t0, const S &t1, S *w ){
	S sum=S(0.0);
	interpolation_cubic_2D_weights( t0, t1, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9], w[10], w[11], w[12], w[13], w[14], w[15] );
	for( int i=0; i<16; i+=4 ){
		sum += v[i+0]*w[i+0];
		sum += v[i+1]*w[i+1];
		sum += v[i+2]*w[i+2];
		sum += v[i+3]*w[i+3];
	}
	return sum;
}

/**
 @brief test function for the 2D cosine interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_cubic_2D(){
	const double t0=0.325f, t1=0.675;
	double o0, o1, o2, o3, v[16], w[16];
	bool result = true;
	
	for( int i=0; i<16; i++ ){
		v[i] = drand48()*5.0;
	}
	
	// check API consistency
	o0 = interpolation_cubic_2D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15], t0, t1 );
	o1 = interpolation_cubic_2D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15], t0, t1, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7], w[8], w[9], w[10], w[11], w[12], w[13], w[14], w[15] );
	o2 = interpolation_cubic_2D( v, t0, t1 );
	o3 = interpolation_cubic_2D( v, t0, t1, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// check input reproduction
	result &= fabs( interpolation_cubic_2D( v, 0.0, 0.0 ) - v[ 5] ) < 1e-9;
	result &= fabs( interpolation_cubic_2D( v, 1.0, 0.0 ) - v[ 6] ) < 1e-9;
	result &= fabs( interpolation_cubic_2D( v, 0.0, 1.0 ) - v[ 9] ) < 1e-9;
	result &= fabs( interpolation_cubic_2D( v, 1.0, 1.0 ) - v[10] ) < 1e-9;

	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D Linear Interpolation (Tensor product formulation) ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for 3D trilinear interpolation
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w000 weight for [x,y,z] = [0,0,0]
 @param[out] w100 weight for [x,y,z] = [1,0,0]
 @param[out] w010 weight for [x,y,z] = [0,1,0]
 @param[out] w110 weight for [x,y,z] = [1,1,0]
 @param[out] w001 weight for [x,y,z] = [0,0,1]
 @param[out] w101 weight for [x,y,z] = [1,0,1]
 @param[out] w011 weight for [x,y,z] = [0,1,1]
 @param[out] w111 weight for [x,y,z] = [1,1,1]
*/
template< typename S >
inline void interpolation_lerp_3D_weights( const S &t0, const S &t1, const S &t2, S &w000, S &w100, S &w010, S &w110, S &w001, S &w101, S &w011, S &w111 ){
	S wlo, whi;
	interpolation_lerp_2D_weights( t0, t1, w000, w100, w010, w110 );
	interpolation_lerp_1D_weights( t2, wlo, whi );
	w001=w000; w101=w100; w011=w010; w111=w110;
	w000*=wlo; w100*=wlo; w010*=wlo; w110*=wlo;
	w001*=whi; w101*=whi; w011*=whi; w111*=whi;
}

/**
 @brief performs trilinear interpolation, using weights computed
 by interpolation_lerp_3D_weights()
 @param[in] v000 value for [x,y,z] = [0,0,0]
 @param[in] v100 value for [x,y,z] = [1,0,0]
 @param[in] v010 value for [x,y,z] = [0,1,0]
 @param[in] v110 value for [x,y,z] = [1,1,0]
 @param[in] v001 value for [x,y,z] = [0,0,1]
 @param[in] v101 value for [x,y,z] = [1,0,1]
 @param[in] v011 value for [x,y,z] = [0,1,1]
 @param[in] v111 value for [x,y,z] = [1,1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_3D( const V &v000, const V &v100, const V &v010, const V &v110, const V &v001, const V &v101, const V &v011, const V &v111, const S &t0, const S &t1, const S &t2 ){
	S w000, w100, w010, w110, w001, w101, w011, w111;
	interpolation_lerp_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v000*w000 + v100*w100 + v010*w010 + v110*w110 
	+ v001*w001 + v101*w101 + v011*w011 + v111*w111;
}

/**
 @brief performs trilinear interpolation, using weights computed
 by interpolation_lerp_3D_weights(), returning the weights
 @param[in] v000 value for [x,y,z] = [0,0,0]
 @param[in] v100 value for [x,y,z] = [1,0,0]
 @param[in] v010 value for [x,y,z] = [0,1,0]
 @param[in] v110 value for [x,y,z] = [1,1,0]
 @param[in] v001 value for [x,y,z] = [0,0,1]
 @param[in] v101 value for [x,y,z] = [1,0,1]
 @param[in] v011 value for [x,y,z] = [0,1,1]
 @param[in] v111 value for [x,y,z] = [1,1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w000 weight for [x,y,z] = [0,0,0]
 @param[out] w100 weight for [x,y,z] = [1,0,0]
 @param[out] w010 weight for [x,y,z] = [0,1,0]
 @param[out] w110 weight for [x,y,z] = [1,1,0]
 @param[out] w001 weight for [x,y,z] = [0,0,1]
 @param[out] w101 weight for [x,y,z] = [1,0,1]
 @param[out] w011 weight for [x,y,z] = [0,1,1]
 @param[out] w111 weight for [x,y,z] = [1,1,1]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_3D( const V &v000, const V &v100, const V &v010, const V &v110, const V &v001, const V &v101, const V &v011, const V &v111, const S &t0, const S &t1, const S &t2, S &w000, S &w100, S &w010, S &w110, S &w001, S &w101, S &w011, S &w111 ){
	interpolation_lerp_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v000*w000 + v100*w100 + v010*w010 + v110*w110 
	+ v001*w001 + v101*w101 + v011*w011 + v111*w111;
}

/**
 @brief performs trilinear interpolation using weights computed
 by interpolation_lerp_3D_weights(), taking input values as an array
 @param[in]  v   array of input value v = [ v000, v100, v010, v110, v001, v101, v011, v111 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_3D( const V *v, const S &t0, const S &t1, const S &t2 ){
	S w000, w100, w010, w110, w001, w101, w011, w111;
	interpolation_lerp_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v[0]*w000 + v[1]*w100 + v[2]*w010 + v[3]*w110 
	+ v[4]*w001 + v[5]*w101 + v[6]*w011 + v[7]*w111;
}

/**
 @brief performs trilinear interpolation using weights computed
 by interpolation_lerp_3D_weights() and returning weights, taking 
 input values and weights as arrays
 @param[in]  v   array of input value v = [ v000, v100, v010, v110, v001, v101, v011, v111 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w   array out output weights in same order as @a v
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_lerp_3D( const V *v, const S &t0, const S &t1, const S &t2, S *w ){
	interpolation_lerp_3D_weights( t0, t1, t2, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7] );
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] + v[3]*w[3] 
	+ v[4]*w[4] + v[5]*w[5] + v[6]*w[6] + v[7]*w[7];
}

/**
 @brief test function for the 2D linear interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_lerp_3D(){
	const double t0=0.325f, t1=0.675, t2=0.156;
	double o0, o1, o2, o3, v[8], w[8];
	bool result = true;
	for( int i=0; i<8; i++ ){
		v[i] = drand48()*5.0;
	}
	
	// test API consistency
	o0 = interpolation_lerp_3D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], t0, t1, t2 );
	o1 = interpolation_lerp_3D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], t0, t1, t2, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7] );
	o2 = interpolation_lerp_3D( v, t0, t1, t2 );
	o3 = interpolation_lerp_3D( v, t0, t1, t2, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// test value reproduction
	result &= fabs( interpolation_lerp_3D( v, 0.0, 0.0, 0.0 ) - v[0] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 1.0, 0.0, 0.0 ) - v[1] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 0.0, 1.0, 0.0 ) - v[2] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 1.0, 1.0, 0.0 ) - v[3] ) < 1e-9;
	
	result &= fabs( interpolation_lerp_3D( v, 0.0, 0.0, 1.0 ) - v[4] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 1.0, 0.0, 1.0 ) - v[5] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 0.0, 1.0, 1.0 ) - v[6] ) < 1e-9;
	result &= fabs( interpolation_lerp_3D( v, 1.0, 1.0, 1.0 ) - v[7] ) < 1e-9;
	
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D Cosine Interpolation (Tensor product formulation) ///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for 3D cosine interpolation
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w000 weight for [x,y,z] = [0,0,0]
 @param[out] w100 weight for [x,y,z] = [1,0,0]
 @param[out] w010 weight for [x,y,z] = [0,1,0]
 @param[out] w110 weight for [x,y,z] = [1,1,0]
 @param[out] w001 weight for [x,y,z] = [0,0,1]
 @param[out] w101 weight for [x,y,z] = [1,0,1]
 @param[out] w011 weight for [x,y,z] = [0,1,1]
 @param[out] w111 weight for [x,y,z] = [1,1,1]
 */
template< typename S >
inline void interpolation_cosine_3D_weights( const S &t0, const S &t1, const S &t2, S &w000, S &w100, S &w010, S &w110, S &w001, S &w101, S &w011, S &w111 ){
	S wlo, whi;
	interpolation_cosine_2D_weights( t0, t1, w000, w100, w010, w110 );
	interpolation_cosine_1D_weights( t2, wlo, whi );
	w001=w000; w101=w100; w011=w010; w111=w110;
	w000*=wlo; w100*=wlo; w010*=wlo; w110*=wlo;
	w001*=whi; w101*=whi; w011*=whi; w111*=whi;
}

/**
 @brief performs 3D cosine interpolation, using weights computed
 by interpolation_lerp_3D_weights()
 @param[in] v000 value for [x,y,z] = [0,0,0]
 @param[in] v100 value for [x,y,z] = [1,0,0]
 @param[in] v010 value for [x,y,z] = [0,1,0]
 @param[in] v110 value for [x,y,z] = [1,1,0]
 @param[in] v001 value for [x,y,z] = [0,0,1]
 @param[in] v101 value for [x,y,z] = [1,0,1]
 @param[in] v011 value for [x,y,z] = [0,1,1]
 @param[in] v111 value for [x,y,z] = [1,1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_3D( const V &v000, const V &v100, const V &v010, const V &v110, const V &v001, const V &v101, const V &v011, const V &v111, const S &t0, const S &t1, const S &t2 ){
	S w000, w100, w010, w110, w001, w101, w011, w111;
	interpolation_cosine_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v000*w000 + v100*w100 + v010*w010 + v110*w110 
	+ v001*w001 + v101*w101 + v011*w011 + v111*w111;
}

/**
 @brief performs 3D cosine interpolation, using weights computed
 by interpolation_cosine_3D_weights(), returning the weights
 @param[in] v000 value for [x,y,z] = [0,0,0]
 @param[in] v100 value for [x,y,z] = [1,0,0]
 @param[in] v010 value for [x,y,z] = [0,1,0]
 @param[in] v110 value for [x,y,z] = [1,1,0]
 @param[in] v001 value for [x,y,z] = [0,0,1]
 @param[in] v101 value for [x,y,z] = [1,0,1]
 @param[in] v011 value for [x,y,z] = [0,1,1]
 @param[in] v111 value for [x,y,z] = [1,1,1]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w000 weight for [x,y,z] = [0,0,0]
 @param[out] w100 weight for [x,y,z] = [1,0,0]
 @param[out] w010 weight for [x,y,z] = [0,1,0]
 @param[out] w110 weight for [x,y,z] = [1,1,0]
 @param[out] w001 weight for [x,y,z] = [0,0,1]
 @param[out] w101 weight for [x,y,z] = [1,0,1]
 @param[out] w011 weight for [x,y,z] = [0,1,1]
 @param[out] w111 weight for [x,y,z] = [1,1,1]
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_3D( const V &v000, const V &v100, const V &v010, const V &v110, const V &v001, const V &v101, const V &v011, const V &v111, const S &t0, const S &t1, const S &t2, S &w000, S &w100, S &w010, S &w110, S &w001, S &w101, S &w011, S &w111 ){
	interpolation_cosine_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v000*w000 + v100*w100 + v010*w010 + v110*w110 
	+ v001*w001 + v101*w101 + v011*w011 + v111*w111;
}

/**
 @brief performs 3D cosine interpolation using weights computed
 by interpolation_cosine_3D_weights(), taking input values as an array
 @param[in]  v   array of input value v = [ v000, v100, v010, v110, v001, v101, v011, v111 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_3D( const V *v, const S &t0, const S &t1, const S &t2 ){
	S w000, w100, w010, w110, w001, w101, w011, w111;
	interpolation_cosine_3D_weights( t0, t1, t2, w000, w100, w010, w110, w001, w101, w011, w111 );
	return v[0]*w000 + v[1]*w100 + v[2]*w010 + v[3]*w110 
	+ v[4]*w001 + v[5]*w101 + v[6]*w011 + v[7]*w111;
}

/**
 @brief performs 3D cosine interpolation using weights computed
 by interpolation_cosine_3D_weights() and returning weights, taking 
 input values and weights as arrays
 @param[in]  v   array of input value v = [ v000, v100, v010, v110, v001, v101, v011, v111 ]
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w   array out output weights in same order as @a v
 @return interpolated value
 */
template< typename V, typename S >
inline V interpolation_cosine_3D( const V *v, const S &t0, const S &t1, const S &t2, S *w ){
	interpolation_cosine_3D_weights( t0, t1, t2, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7] );
	return v[0]*w[0] + v[1]*w[1] + v[2]*w[2] + v[3]*w[3] 
	+ v[4]*w[4] + v[5]*w[5] + v[6]*w[6] + v[7]*w[7];
}

/**
 @brief test function for the 3D cosine interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_cosine_3D(){
	const double t0=0.325f, t1=0.675, t2=0.156;
	double o0, o1, o2, o3, v[8], w[8];
	bool result = true;
	for( int i=0; i<8; i++ ){
		v[i] = drand48()*5.0;
	}
	o0 = interpolation_cosine_3D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], t0, t1, t2 );
	o1 = interpolation_cosine_3D( v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], t0, t1, t2, w[0], w[1], w[2], w[3], w[4], w[5], w[6], w[7] );
	o2 = interpolation_cosine_3D( v, t0, t1, t2 );
	o3 = interpolation_cosine_3D( v, t0, t1, t2, w );
	result &= fabs(o0-o1) < 1e-9 && fabs(o0-o2) < 1e-9 && fabs(o0-o3) < 1e-9;
	
	// test value reproduction
	result &= fabs( interpolation_cosine_3D( v, 0.0, 0.0, 0.0 ) - v[0] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 1.0, 0.0, 0.0 ) - v[1] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 0.0, 1.0, 0.0 ) - v[2] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 1.0, 1.0, 0.0 ) - v[3] ) < 1e-9;
	
	result &= fabs( interpolation_cosine_3D( v, 0.0, 0.0, 1.0 ) - v[4] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 1.0, 0.0, 1.0 ) - v[5] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 0.0, 1.0, 1.0 ) - v[6] ) < 1e-9;
	result &= fabs( interpolation_cosine_3D( v, 1.0, 1.0, 1.0 ) - v[7] ) < 1e-9;
	
	return result;	
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// 3D Cubic Interpolation (Tensor product formulation) ////////////////////////////////////////////////
// This uses a slighyly different API than the others, due to the large numbers of points /////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 @brief computes the weights needed for cubic interpolation in 3D
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @param[out] w   array of weights ordered lexographically by first, then second, then third axis
*/
template< typename S >
inline void interpolation_cubic_3D_weights( const S &t0, const S &t1, const S &t2, S *w ){
	S wx[4], wy[4], wz[4];
	interpolation_cubic_1D_weights( t0, wx[0], wx[1], wx[2], wx[3] );
	interpolation_cubic_1D_weights( t1, wy[0], wy[1], wy[2], wy[3] );
	interpolation_cubic_1D_weights( t2, wz[0], wz[1], wz[2], wz[3] );
	int wid = 0;
	// loop order is important here....third axis must be outer loop
	// first axis must be inner (unrolled) loop.
	for( int k=0; k<4; k++ ){
		for( int j=0; j<4; j++ ){
			// unroll the loop for a modicum of efficiency
			w[wid++] = wx[0]*wy[j]*wz[k];
			w[wid++] = wx[1]*wy[j]*wz[k];
			w[wid++] = wx[2]*wy[j]*wz[k];
			w[wid++] = wx[3]*wy[j]*wz[k];
		}
	}
}

/**
 @brief performs 3D cubic interpolation using weights computed
 by interpolation_cubic_3D_weights(), taking input values as an array
 @param[in]  v   array of input values, lexographically ordered by first then second then third axis
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value 
*/
template< typename V, typename S >
inline V interpolation_cubic_3D( const V *v, const S &t0, const S &t1, const S &t2 ){
	S sum=S(0.0), w[64];
	interpolation_cubic_3D_weights( t0, t1, t2, w );
	for( int i=0; i<64; i+= 8 ){
		sum += v[i+0]*w[i+0];
		sum += v[i+1]*w[i+1];
		sum += v[i+2]*w[i+2];
		sum += v[i+3]*w[i+3];
		sum += v[i+4]*w[i+4];
		sum += v[i+5]*w[i+5];
		sum += v[i+6]*w[i+6];
		sum += v[i+7]*w[i+7];
	}
	return sum;
}

/**
 @brief performs 3D cubic interpolation using weights computed
 by interpolation_cubic_3D_weights(), taking input values as an array
 @param[in]  v   array of input values, lexographically ordered by first then second then third axis
 @param[in]  t0  fractional distance along first axis
 @param[in]  t1  fractional distance along second axis
 @param[in]  t2  fractional distance along third axis
 @return interpolated value 
 */
template< typename V, typename S >
inline V interpolation_cubic_3D( const V *v, const S &t0, const S &t1, const S &t2, S *w ){
	S sum=S(0.0);
	interpolation_cubic_3D_weights( t0, t1, t2, w );
	for( int i=0; i<64; i+= 8 ){
		sum += v[i+0]*w[i+0];
		sum += v[i+1]*w[i+1];
		sum += v[i+2]*w[i+2];
		sum += v[i+3]*w[i+3];
		sum += v[i+4]*w[i+4];
		sum += v[i+5]*w[i+5];
		sum += v[i+6]*w[i+6];
		sum += v[i+7]*w[i+7];
	}
	return sum;
}

/**
 @brief test function for the 3D cubic interpolation routines, currently just
 checks that the methods are consistent, i.e. that they all produce the same
 output for the same input
 @return true if test passes false otherwise
 */
static bool interpolation_test_cubic_3D(){
	const double t0=0.325f, t1=0.675, t2=0.156;
	double o0, o1, v[64], w[64];
	bool result = true;
	for( int i=0; i<64; i++ ){
		v[i] = drand48()*5.0;
	}
	
	// test API consistency
	o0 = interpolation_cubic_3D( v, t0, t1, t2 );
	o1 = interpolation_cubic_3D( v, t0, t1, t2, w );
	result &= fabs(o0-o1) < 1e-9;
	
	// test value reproduction
	result &= fabs( interpolation_cubic_3D( v, 0.0, 0.0, 0.0 ) - v[21] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 1.0, 0.0, 0.0 ) - v[22] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 0.0, 1.0, 0.0 ) - v[25] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 1.0, 1.0, 0.0 ) - v[26] ) < 1e-9;
	
	result &= fabs( interpolation_cubic_3D( v, 0.0, 0.0, 1.0 ) - v[37] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 1.0, 0.0, 1.0 ) - v[38] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 0.0, 1.0, 1.0 ) - v[41] ) < 1e-9;
	result &= fabs( interpolation_cubic_3D( v, 1.0, 1.0, 1.0 ) - v[42] ) < 1e-9;
	
	return result;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
// Test driver ////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
static bool interpolation_test(){
	bool result = true;
	// 1D tests
	result &= interpolation_test_lerp_1D();
	result &= interpolation_test_cosine_1D();
	result &= interpolation_test_cubic_1D();
	
	// 2D tests
	result &= interpolation_test_lerp_2D();
	result &= interpolation_test_cosine_2D();
	result &= interpolation_test_cubic_2D();
	
	// 3D tests
	result &= interpolation_test_lerp_3D();
	result &= interpolation_test_cosine_2D();
	result &= interpolation_test_cubic_3D();
	
	return result;
}

#endif
