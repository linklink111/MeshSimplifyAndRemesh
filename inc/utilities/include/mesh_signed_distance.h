#ifndef MESH_SIGNED_DISTANCE_H
#define MESH_SIGNED_DISTANCE_H

/**
 @file mesh_signed_distance.h
 @author James Gregson (james.gregson@gmail.com)
 @copyright Copyright (C) 2014 James Gregson
 
 Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @brief Generic function to compute narrow-band signed distance fields of triangle meshes based on regular grids based on the angle-weighted normals of Bærentzen and Aanæs, 2002.  The function is fully generic and depends only on the C++ standard library as well as a user-provided 3-component vector and 3D image classes that must conform to a minimally intrusive interface.
 
	The vector class must provide the following methods (a,b are vectors, s, x, y, z are scalars):
		- default constructor: a() initializing components to zero
		- constructor: a(x,y,z) initializing components to x,y,z
		- indexing operator: a[i], i=0,1,2 returning writable reference to the i'th component
		- arithmetic operations: a+b, a-b, a+=b, a-=b for vector addition/subtraction
		- scaling operations: a*s, b*s, a*=s, b*=s for multiplication/division by a scalar
		
	These operations are commonly found in 3-component vector classes so it is expected that these requirements are not difficult to meet in practice.

	The image class must provide the call operator img(i,j,k) that returns a writeable reference to image components, i.e. img(i,j,k) must return a writable reference to the grid point [i,j,k].  The spatial dimensions of the image to be calculated are specified by an axis aligned bounding box AABB packed as a 6-tuple of real values [xmin, xmax, ymin, ymax, zmin, zmax] and the points along each dimension are assumed to be uniformly spaced such the first and last point along an axis correspond to the min and max values in the AABB, i.e. img(0,0,0) correponds to a point with spatial location [xmin,ymin,zmin] and img(nx-1, ny-1, nz-1) corresponds to point [xmax,ymax,zmax].
	
	To generate a distance field, simply call the function:
	
	template< typename image, typename real3, typename real, typename index=int >
	void mesh_to_signed_distance( const std::vector<real3>& vtx, const std::vector<index>& tri, const real* aabb, const index* dim, image& signed_dis, const real band_size );
	
	with two std::vectors, one containing mesh vertices and the other containing flat-packed triangle vertex indices (zero-based).  The aabb array is a 6-tuple packed as [xmin,xmax,ymin,ymax,zmin,zmax], signed_dist is an image to hold the signed distance values.  dim is an 3-tuple [nx,ny,nz] specifying the number of points in each direction while band-size sets the maximum absolute distance from the mesh to compute signed distance values.
*/

#include<map>
#include<vector>
#include<cmath>

namespace mesh_signed_distance {
	
	/**
	 @brief Function used in error checking, prints/throws the error if condition evaluates to false. Used for input checking.  checks are minimally intrusive but can be disabled entirely by defining the preprocessor macro MESH_SIGNED_DISTANCE_DISABLE_CHECKS
	*/
	inline void check( const bool condition, const char *msg ){
		if( !condition ){
			std::cout << msg << std::endl;
			throw msg;
		}
	}
	
	/**
	 @brief returns a real3 with the component-wise maximum of a and b
	*/
	template< typename real3 >
	inline real3 max( const real3& a, const real3& b ){
		return real3( std::max(a[0],b[0]), std::max(a[1],b[1]), std::max(a[2],b[2]) );
	}
	
	/**
	 @brief returns a real3 with the component-wise minimum of a and b
	*/
	template< typename real3 >
	inline real3 min( const real3& a, const real3& b ){
		return real3( std::min(a[0],b[0]), std::min(a[1],b[1]), std::min(a[2],b[2]) );
	}
	
	/**
	 @brief returns the dot-product of a and b
	*/
	template< typename real3, typename real=double >
	inline real dot( const real3& a, const real3& b ){
		return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
	}
	
	/**
	 @brief returns the cross-product of a and b
	*/
	template< typename real3 >
	inline real3 cross( const real3& a, const real3& b ){
		return real3( a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0] );
	}
	
	/**
	 @brief returns the squared length of a
	*/
	template< typename real3, typename real=double >
	inline real length_squared( const real3& a ){
		return dot(a,a);
	}
	
	/**
	 @brief returns the lenght of a
	*/
	template< typename real3, typename real=double >
	inline real length( const real3& a ){
		return (real)sqrt( length_squared<real3,real>(a) );
	}
	
	/**
	 @brief safely normalizes the vector a in place, returning the length. if the length of a is less than the zero_thresh value, uses zero_thresh in place of the length to prevent divide-by-zero
	*/
	template< typename real3, typename real=double >
	inline real safe_normalize( real3& a, const real zero_thresh=real(1e-5) ){
		real l = std::max( length(a), zero_thresh );
		a /= l;
		return l;
	}
	
	/**
	 @brief compute the internal angles of a triangle, based on the vertex positions A, B, C. results are packed in order, e.g. return_value[0] is the internal angle at point A, return_value[1] is internal angle at B
	*/
	template< typename real3, typename real >
	inline real3 triangle_angles( const real3& A, const real3& B, const real3& C ){
		real a = length<real3,real>(C-B);
		real b = length<real3,real>(C-A);
		real c = length<real3,real>(B-A);
		real alpha = acos( (b*b + c*c - a*a)/(2.0*b*c) );
		real beta  = acos( (a*a + c*c - b*b)/(2.0*a*c) );
		real gamma = M_PI - alpha - beta; //acos( (a*a + b*b - c*c)/(2.0*a*b) );
		return real3( alpha, beta, gamma );
	}
	
	/**
	 @brief returns the axis-aligned bounding box of the triangle ABC, storing the near lower left corner in minim and the far upper right corner in maxim
	*/
	template< typename real3 >
	inline void triangle_bbox( const real3& A, const real3& B, const real3& C, real3& minim, real3& maxim ){
		minim = maxim = A;
		minim = min( minim, min(B,C) );
		maxim = max( maxim, max(B,C) );
	}
	
	/**
	 @brief transforms points from world coordinate to grid coordinates, with no range checking performed. aabb is the axis aligned bounding box for the domain stored as [xmin,xmax,ymin,ymax,zmin,zmax], dim is the number of points in the x, y and z directions. the first and last points are map to the AABB min and max values for the axis.
	*/
	template< typename real3, typename real, typename index >
	inline real3 _world_to_grid( const real* aabb, const index* dim, const real3& w ){
		return real3( real(dim[0]-1)*(w[0]-aabb[0])/(aabb[1]-aabb[0]), real(dim[1]-1)*(w[1]-aabb[2])/(aabb[3]-aabb[2]), real(dim[2]-1)*(w[2]-aabb[4])/(aabb[5]-aabb[4]) );
	}
	
	/**
	 @brief performs the inverse operation of _world_to_grid, tranforming points from grid coordinates to world coordinates.
	*/
	template< typename real3, typename real, typename index >
	inline real3 _grid_to_world( const real* aabb, const index* dim, const real3& g ){
		return real3( g[0]*(aabb[1]-aabb[0])/real(dim[0]-1)+aabb[0], g[1]*(aabb[3]-aabb[2])/real(dim[1]-1)+aabb[2], g[2]*(aabb[5]-aabb[4])/real(dim[2]-1)+aabb[4] );
	}
	
	/**
	 @brief enumeration of the seven possible triangle features that a point can be closest to, there are three vertices, three edges and one facce.
	*/
	enum {
		TRIANGLE_VERTEX_0,
		TRIANGLE_VERTEX_1,
		TRIANGLE_VERTEX_2,
		TRIANGLE_EDGE_01,
		TRIANGLE_EDGE_12,
		TRIANGLE_EDGE_20,
		TRIANGLE_FACE_012,
	};
	
	/**
	 @brief helper class to store the return type of distance queries against a triangle mesh, each instance encodes the closest point and the type of the point (one of the seven possible closest features) in order to allow the appropriate angle-weighted normal value to be looked up and used for sidedness testing.
	*/
	template< typename real3 >
	class triangle_point_type {
	private:
		real3		m_point;
		int			m_type;
	public:
		/** default constructor, do nothing */
		triangle_point_type(){
		}
		/** constuct based on the closest point and the type */
		triangle_point_type( const real3& point, const int& type ) : m_point(point), m_type(type) {
		}
		
		/** cast to real3 to allow seamless arithmetic operations */
		inline operator const real3&() const {
			return m_point;
		}
		/** get the closest point */
		inline const real3& point() const {
			return m_point;
		}
		/** get the type of the closest point */
		inline const int& type() const {
			return m_type;
		}
	};
	
	/**
	 @brief function to return the closest point on a line-segment, returned as a triangle_point_type. tags for the three possible cases (start vertex, end vertex, edge interior) are passed in as arguments and the appropriate tag returned to allow selection of the proper angle-weighted normal value
	*/
	template< typename real3, typename real=double >
	inline triangle_point_type<real3> closest_point_on_segment( const real3& S, const real3& E, const real3& p, const int s_type=TRIANGLE_VERTEX_0, const int e_type=TRIANGLE_VERTEX_1, const int edge_type=TRIANGLE_EDGE_01 ){
		real3 n = E-S;
		real3 d = p-S;
		real t = dot(d,n)/dot(n,n); //d.dot(n)/n.dot(n);
		if( t <= 0.0 ){
			return triangle_point_type<real3>( S, s_type );
		} else if( t >= 1.0 ){
			return triangle_point_type<real3>( E, e_type );
		} else {
			return triangle_point_type<real3>( S + n*t, edge_type );
		}
	}
	
	/**
	 @brief computes the barycentric coordinates of the closest point in the plane of triangle T1,T2,T3 to the test point p using a least-squares approach. det_thresh is used to detect degenerate triangles (which only occur when one of the three edges is zero length or two edges are parallel), in which case the system is modified slightly by the addition of a small identity term.
	*/
	template< typename real3, typename real=double >
	inline real3 triangle_barycentric_coordinates( const real3& T0, const real3& T1, const real3 &T2, const real3& p, const real det_thresh=real(1e-5) ){
		// compute two edge vectors and a delta vector,
		// transforming the triangle to have its third
		// vertex at the origin
		const real3 e0 = T0 - T2;
		const real3 e1 = T1 - T2;
		const real3 ep = p  - T2;
		
		// build a least-squares system for the barycentric
		// first two barycentric coordinates minimizing the
		// squared distance to the point. the third coordinate
		// is eliminated from the system by the subtraction
		// of T2 when forming the edge and delta vectors
		real AtA[][2] = {
			{ dot(e0,e0), dot(e0,e1) },
			{ dot(e1,e0), dot(e1,e1) }
		};
		real b[] = { dot(e0,ep), dot(e1,ep) };
		
		// solve the system by Crout's method, first compute
		// the determinant of the matrix AtA, then multiply
		// the right hand-side by the inverse of AtA,
		// computing the product on the fly rather than
		// explicitly forming AtA^(-1)
		real det = AtA[0][0]*AtA[1][1] - AtA[0][1]*AtA[1][0];
		if( fabs(det) < det_thresh ){
			// detect degenerate triangles and modify the the least squares
			// system by adding a small factor of identity, i.e. use
			// Tihkonov regularization to try to get a reasonable solution.
			// then recalculate the determinant
			AtA[0][0] += 1e-4;
			AtA[1][1] += 1e-4;
			det = AtA[0][0]*AtA[1][1] - AtA[0][1]*AtA[1][0];
		}
		// now multiply through by the inverse
		const real x[] = {
			(+AtA[1][1]*b[0] - AtA[0][1]*b[1])/det,
			(-AtA[1][0]*b[0] + AtA[0][0]*b[1])/det
		};
		
		// return the barycentric coordinates, computing
		// the third component from the condition that
		// x[0]+x[1]+x[2] = 1
		return real3( x[0], x[1], 1.0-x[0]-x[1] );
	}
	
	/**
	 @brief returns the closest point and closest point type on a triangle using a case-based method.  First barycentric coordinates are computed for the closest point in the plane of the triangle to the test point. These are then tested to see if result is interior to the triangle (which is returned if true). Otherwise the edges are tested based on the barycentric coordinates and the minimum distance result (and type) returned.  det_thresh is a factor determining when a triangle is degenerate based on the determinant of the least-squares system used to compute the barycentric coordinate of the closest point in the plane of the triangle while bary_thresh is a threshold used to distinguish between edge and interior intersections.
	*/
	template< typename real3, typename real=double >
	inline triangle_point_type<real3> closest_point_on_triangle( const real3& T0, const real3& T1, const real3& T2, const real3& p, const real det_thresh=real(1e-5), const real bary_thresh=real(1e-4) ){
		// compute the barycentric coordinates of the in-plane point closest to the test poin p
		real3 bary = triangle_barycentric_coordinates( T0, T1, T2, p, det_thresh );
		if( bary[0] > bary_thresh && bary[1] > bary_thresh && bary[2] > bary_thresh && fabs(real(1.0)-bary[0]-bary[1]-bary[2]) < bary_thresh ){
			return triangle_point_type<real3>( T0*bary[0] + T1*bary[1] + T2*bary[2], TRIANGLE_FACE_012 );
		}
		
		// point does not map to interior of face, check the barycentric
		// coordinate to see which edges the point is on the incorrect side
		// of and compute the closest point on all violated edges
		real						test_dis, best_dis = 1e10;
		triangle_point_type<real3>	test_pnt, best_pnt;
		if( bary[0] <= bary_thresh ){
			test_pnt = closest_point_on_segment<real3,real>( T1, T2, p, TRIANGLE_VERTEX_1, TRIANGLE_VERTEX_2, TRIANGLE_EDGE_12 );
			if( (test_dis=length_squared(p-test_pnt)) < best_dis ){
				best_dis = test_dis;
				best_pnt = test_pnt;
			}
		}
		if( bary[1] <= bary_thresh ){
			test_pnt = closest_point_on_segment<real3,real>( T2, T0, p, TRIANGLE_VERTEX_2, TRIANGLE_VERTEX_0, TRIANGLE_EDGE_20 );
			if( (test_dis=length_squared(p-test_pnt)) < best_dis ){
				best_dis = test_dis;
				best_pnt = test_pnt;
			}
		}
		if( bary[2] <= bary_thresh ){
			test_pnt = closest_point_on_segment<real3,real>( T0, T1, p, TRIANGLE_VERTEX_0, TRIANGLE_VERTEX_1, TRIANGLE_EDGE_01 );
			if( (test_dis=length_squared(p-test_pnt)) < best_dis ){
				best_dis = test_dis;
				best_pnt = test_pnt;
			}
		}
		// return the best of the found points
		return best_pnt;
	}
	
	/**
	 @brief sole public function of the library, takes a mesh structure passed in as a vector of vertex positions and a vector of flat-packed 0-based triangle indices and computes the signed distance to the mesh in a narrow band of size band_size on the image signed_dis.  Template parameters real3 are 3-component vectors which must supply a basic interface described at the top of this header.  The aabb and dim arrays specific the bounding box (packed as a 6-tuple [xmin,xmax,ymin,ymax,zmin,zmax]) and dimensions (a 3-tuple [nx,ny,nz]).  The first and last points along an axis correspond to the min and max values from the AABB entries for that axis.
	 
		The method operates by computing angle-weighted normals for the triangles, vertices and edges then looping over the triangles and, for all image vertices within the band_size units of the triangle bounding box, computing the minimum distance to the triangle.  Thus the method will be most efficient when i) triangles are well-shaped and axis aligned, ii) the image resolution is comparable to the surface mesh resolution and iii) the band-size is kept small so relatively few points are within band_size units of the surface.  Image points outside of the narrow band are untouched and can be initialized to +inf in order to allow the output of this function to serve as the seed points for a fast-marching or fast-sweeping method to complete the domain.
		
		This could be improved by determining the best 2D projection of each triangle and then only looping over image vertices that are within the band-size of the triangle in the remaining direction, to avoid traversing many vertices outside the band for large, off-axis triangles. However 256^3 images can be generated for ~70k triangle models in just a few seconds currently so this is not an urgent update
	*/
	template< typename image, typename real3, typename real, typename index=int >
	void mesh_to_signed_distance( const std::vector<real3>& vtx, const std::vector<index>& tri, const real* aabb, const index* dim, image& signed_dis, const real band_size ){
		// begin by initializing arrays of triangle, vertex and edge normals
		// using the angle-weighted normals of Bærentzen and Aanæs, 2002
		typedef std::pair<int,int>	iipair;
		std::vector<real3>			nrm_vtx( vtx.size() );
		std::vector<real3>			nrm_tri( tri.size()/3 );
		std::map< iipair, real3 >	nrm_edg;  // TODO: replace with unordered_map, but needs hash function defined for iipair
		
#ifndef MESH_SIGNED_DISTANCE_DISABLE_CHECKS
		// check that the number of entries in the triangles list is a multiple of 3
		check( tri.size()%3 == 0, "Error, triangle index array not a multiple of 3" );

		// check that grid points as computed by the _grid_to_world function agree with the bounding
		// box extents for both the near lower left and far upper right corners
		check( length_squared(_grid_to_world(aabb,dim,real3(0,0,0))-real3(aabb[0],aabb[2],aabb[4])) < 1e-4, "Error in _grid_to_world()" );
		check( length_squared(_grid_to_world(aabb,dim,real3(dim[0]-1,dim[1]-1,dim[2]-1))-real3(aabb[1],aabb[3],aabb[5])) < 1e-4, "Error in _grid_to_world()" );
		
		// now check that the function _world_to_grid inverts the operation performed by the function
		// _grid_to_world properly as a sanity check
		check( length_squared(_world_to_grid(aabb,dim,_grid_to_world(aabb,dim,real3(0,0,0)))-real3(0,0,0)) < 1e-4, "Error in _world_to_grid()" );
		check( length_squared(_world_to_grid(aabb,dim,_grid_to_world(aabb,dim,real3(dim[0]-1,dim[1]-1,dim[2]-1)))-real3(dim[0]-1,dim[1]-1,dim[2]-1)) < 1e-4, "Error in _world_to_grid()" );
#endif
		
		// loop over the triangle and compute their normals and angles
		for( int i=0; i<tri.size(); i+=3 ){
			// grab the vertices
			index vid[] = { tri[i+0], tri[i+1], tri[i+2] };
			
#ifndef MESH_SIGNED_DISTANCE_DISABLE_CHECKS
			// check that the vertex ids are valid under zero-based indexing from the triangle array
			check( vid[0] >= 0 && vid[0] < tri.size() && vid[1] >= 0 && vid[1] < tri.size() && vid[2] >= 0 && vid[2] < tri.size(), "Error in triangle vertex indices" );
#endif
			// look up the vertices
			real3 v[] = { vtx[vid[0]], vtx[vid[1]], vtx[vid[2]] };
			
			// grab the normal
			real3 n = cross( v[1]-v[0], v[2]-v[0]);
			safe_normalize( n, 1e-5 );
			
			// update the triangle, edge and vertex normals. triangle and edge
			// normals are simply summed, but the vertex normals are weighted by
			// the incident angle of each triangle to the vertex in question
			real3 angles = triangle_angles<real3,real>( v[0], v[1], v[2] );
			nrm_tri[i/3] = n;
			nrm_edg[ iipair( std::min(vid[0],vid[1]), std::max(vid[0],vid[1]) ) ] += n;
			nrm_edg[ iipair( std::min(vid[1],vid[2]), std::max(vid[1],vid[2]) ) ] += n;
			nrm_edg[ iipair( std::min(vid[2],vid[0]), std::max(vid[2],vid[0]) ) ] += n;
			nrm_vtx[ vid[0] ] += n*angles[0];
			nrm_vtx[ vid[1] ] += n*angles[1];
			nrm_vtx[ vid[2] ] += n*angles[2];
		}
		
		// now perform a second loop over the triangles to set the values
		triangle_point_type<real3> cp;
		const real3 off( band_size, band_size, band_size );
		//const real3 h( (aabb[1]-aabb[0])/real(dim[0]), (aabb[3]-aabb[2])/real(dim[1]), (aabb[5]-aabb[4])/real(dim[2]) );
		for( int i=0; i<tri.size(); i+=3 ){
			real3 minim, maxim, p, n, d;
			const index vid[] = { tri[i+0], tri[i+1], tri[i+2] };
			const real3 v[] = { vtx[vid[0]], vtx[vid[1]], vtx[vid[2]] };
			
			// grab the triangle bounding box
			triangle_bbox( v[0], v[1], v[2], minim, maxim );
			minim = _world_to_grid( aabb, dim, minim-off );
			maxim = _world_to_grid( aabb, dim, maxim+off );
			
			// get the bounding indices in the grid coordinates
			index rbnd[] = { std::max( index(minim[0]), 0 ), std::min(index(maxim[0]+1.0), dim[0]-1) };
			index sbnd[] = { std::max( index(minim[1]), 0 ), std::min(index(maxim[1]+1.0), dim[1]-1) };
			index tbnd[] = { std::max( index(minim[2]), 0 ), std::min(index(maxim[2]+1.0), dim[2]-1) };
			
			// loop over the range of indices to update the grid points
			for( index t=tbnd[0]; t<=tbnd[1]; t++ ){
				for( index s=sbnd[0]; s<=sbnd[1]; s++ ){
					for( index r=rbnd[0]; r<=rbnd[1]; r++ ){
						p = _grid_to_world( aabb, dim, real3(r,s,t) );
						// find the closest signed distance on the triangle
						cp = closest_point_on_triangle( v[0], v[1], v[2], p );
						
						// compute the vector from the closest point to the surface
						// and check if it has a smaller absolute signed distance
						// than the existing value, if so compute the sign and update
						// the stored signed distance value
						d = p-cp;
						if( length_squared(d) < signed_dis(r,s,t)*signed_dis(r,s,t) ){
							// check the 7 cases for the featue that the closest
							// point is located on and choose the appropriate
							// pseudo-normal vector
							switch( cp.type() ){
								case TRIANGLE_FACE_012:
									n = nrm_tri[i/3];
									break;
								case TRIANGLE_EDGE_01:
									n = nrm_edg[ iipair( std::min(vid[0],vid[1]), std::max(vid[0],vid[1]) ) ];
									break;
								case TRIANGLE_EDGE_12:
									n = nrm_edg[ iipair( std::min(vid[1],vid[2]), std::max(vid[1],vid[2]) ) ];
									break;
								case TRIANGLE_EDGE_20:
									n = nrm_edg[ iipair( std::min(vid[2],vid[0]), std::max(vid[2],vid[0]) ) ];
									break;
								case TRIANGLE_VERTEX_0:
									n = nrm_vtx[ vid[0] ];
									break;
								case TRIANGLE_VERTEX_1:
									n = nrm_vtx[ vid[1] ];
									break;
								case TRIANGLE_VERTEX_2:
									n = nrm_vtx[ vid[2] ];
									break;
							}
							
							// update the signed distance value and sign
							signed_dis(r,s,t) = length(d)*( dot(d,n) >= 0.0 ? 1.0 : -1.0 );
						}
					}
				}
			}
		}
	}
	
};
	
#endif