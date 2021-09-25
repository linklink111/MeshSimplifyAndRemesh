#ifndef AABB_H
#define AABB_H

#include<algorithm>
#include<cassert>

#include"vec3.h"

namespace geom {
    
    template< typename real >
    class aabb {
    private:
        real    m_bounds[6];
    public:
        enum {
            XMIN,
            XMAX,
            YMIN,
            YMAX,
            ZMIN,
            ZMAX,
        };
        
        aabb(){
            set_bounds( 0, 0, 0, 0, 0, 0 );
        }
        
        aabb( const real *bounds ){
            set_bounds(bounds);
        }
        
        aabb( const real xmin, const real xmax, const real ymin, const real ymax, const real zmin, const real zmax ){
            set_bounds(xmin,xmax,ymin,ymax,zmin,zmax);
        }
        
        aabb( const real *minim, const real *maxim ){
            set_bounds( minim, maxim );
        }
        
        aabb( const vec3<real> &minim, const vec3<real> &maxim ){
            set_bounds(minim,maxim);
        }
        
        void set_bounds( const real *bounds ){
            set_bounds( bounds[XMIN], bounds[XMAX], bounds[YMIN], bounds[YMAX], bounds[ZMIN], bounds[ZMAX] );
        }
        
        void set_bounds( const real xmin, const real xmax, const real ymin, const real ymax, const real zmin, const real zmax ){
            assert( xmax >= xmin && ymax >= ymin && zmax >= zmin );
            m_bounds[XMIN] = xmin; m_bounds[XMAX] = xmax;
            m_bounds[YMIN] = ymin; m_bounds[YMAX] = ymax;
            m_bounds[ZMIN] = zmin; m_bounds[ZMAX] = zmax;
        }
        
        void set_bounds( const real *minim, const real *maxim ){
            set_bounds( minim[0], maxim[0], minim[1], maxim[1], minim[2], maxim[2] );
        }
        
        void set_bounds( vec3<real> &minim, vec3<real> &maxim ){
            set_bounds( minim.x(),maxim.x(),minim.y(),maxim.y(),minim.z(),maxim.z() );
        }
        
        void get_bounds( real *bounds ) const {
            bounds[XMIN] = m_bounds[XMIN]; bounds[XMAX] = m_bounds[XMAX];
            bounds[YMIN] = m_bounds[YMIN]; bounds[YMAX] = m_bounds[YMAX];
            bounds[ZMIN] = m_bounds[ZMIN]; bounds[ZMAX] = m_bounds[ZMAX];
        }
        
        void get_bounds( real &xmin, real &xmax, real &ymin, real &ymax, real &zmin, real &zmax ) const {
            xmin = m_bounds[XMIN]; xmax = m_bounds[XMAX];
            ymin = m_bounds[YMIN]; ymax = m_bounds[YMAX];
            zmin = m_bounds[ZMIN]; zmax = m_bounds[ZMAX];
        }
        
        void get_bounds( real *minim, real *maxim ) const {
            minim[0] = m_bounds[XMIN]; maxim[0] = m_bounds[XMAX];
            minim[1] = m_bounds[YMIN]; maxim[1] = m_bounds[YMAX];
            minim[2] = m_bounds[ZMIN]; maxim[2] = m_bounds[ZMAX];
        }
        
        void get_bounds( vec3<real> &minim, vec3<real> &maxim ) const {
            minim.x() = m_bounds[XMIN]; maxim.x() = m_bounds[XMAX];
            minim.y() = m_bounds[YMIN]; maxim.y() = m_bounds[YMAX];
            minim.z() = m_bounds[ZMIN]; maxim.z() = m_bounds[ZMAX];
        }
        
        inline real & operator[] ( const int id ) {
            assert( id >= 0 && id < 6 );
            return m_bounds[id];
        }
        
        inline const real & operator[] ( const int id ) const {
            assert( id >= 0 && id < 6 );
            return m_bounds[id];
        }
        
        inline real xmin() const {
            return m_bounds[XMIN];
        }
        
        inline real xmax() const {
            return m_bounds[XMAX];
        }
        
        inline real ymin() const {
            return m_bounds[YMIN];
        }
        
        inline real ymax() const {
            return m_bounds[YMAX];
        }

        inline real zmin() const {
            return m_bounds[ZMIN];
        }
        
        inline real zmax() const {
            return m_bounds[ZMAX];
        }
        
        inline real dx() const {
            return m_bounds[XMAX] - m_bounds[XMIN];
        }
        
        inline real dy() const {
            return m_bounds[YMAX] - m_bounds[YMIN];
        }
        
        inline real dz() const {
            return m_bounds[ZMAX] - m_bounds[ZMIN];
        }
        
        inline real volume() const {
            return dx()*dy()*dz();
        }
        
        inline real surface_area() const {
            const real d[] = { dx(), dy(), dz() };
            return 2*( d[0]*d[1] + d[0]*d[2] + d[1]*d[2] );
        }
        
        inline real vertices( real *verts ) const {
            verts[ 0] = m_bounds[XMIN]; verts[ 1] = m_bounds[YMIN]; verts[ 2] = m_bounds[ZMIN];
            verts[ 3] = m_bounds[XMAX]; verts[ 4] = m_bounds[YMIN]; verts[ 5] = m_bounds[ZMIN];
            verts[ 6] = m_bounds[XMAX]; verts[ 7] = m_bounds[YMIN]; verts[ 8] = m_bounds[ZMAX];
            verts[ 9] = m_bounds[XMIN]; verts[10] = m_bounds[YMIN]; verts[11] = m_bounds[ZMAX];
            verts[12] = m_bounds[XMIN]; verts[13] = m_bounds[YMAX]; verts[14] = m_bounds[ZMIN];
            verts[15] = m_bounds[XMAX]; verts[16] = m_bounds[YMAX]; verts[17] = m_bounds[ZMIN];
            verts[18] = m_bounds[XMAX]; verts[19] = m_bounds[YMAX]; verts[20] = m_bounds[ZMAX];
            verts[21] = m_bounds[XMIN]; verts[22] = m_bounds[YMAX]; verts[23] = m_bounds[ZMAX];
        }
        
        inline bool parametric_coordinates( const real *pnt, real *param ) const {
            param[0] = (pnt[0]-m_bounds[XMIN])/(m_bounds[XMAX]-m_bounds[XMIN]);
            param[1] = (pnt[1]-m_bounds[YMIN])/(m_bounds[YMAX]-m_bounds[YMIN]);
            param[2] = (pnt[2]-m_bounds[ZMIN])/(m_bounds[ZMAX]-m_bounds[ZMIN]);
            return param[0] >= 0.0 && param[0] <= 1.0 && param[1] >= 0.0 && param[1] <= 1.0 && param[2] >= 0.0 && param[2] <= 1.0;
        }
        
        inline bool intersects_point( const real *pnt ) const {
            return pnt[0] >= m_bounds[XMIN] && pnt[0] <= m_bounds[XMAX] && pnt[1] >= m_bounds[YMIN] && pnt[1] <= m_bounds[YMAX] && pnt[2] >= m_bounds[ZMIN] && pnt[2] <= m_bounds[ZMAX];
        }
        
        inline bool intersects_sphere( const real *center, const real radius ) const {
            real c[3], d[3];
            c[0] = std::max( m_bounds[XMIN], std::min( m_bounds[XMAX], center[0] ) );
            c[1] = std::max( m_bounds[YMIN], std::min( m_bounds[YMAX], center[1] ) );
            c[2] = std::max( m_bounds[ZMIN], std::min( m_bounds[ZMAX], center[2] ) );
            d[0] = center[0]-c[0];
            d[1] = center[1]-c[1];
            d[2] = center[2]-c[2];
            return d[0]*d[0] + d[1]*d[1] + d[2]*d[2] <= radius*radius;
        }
        
        inline bool intersects_aabb( const aabb &box ) const {
            return m_bounds[XMIN] <= box.m_bounds[XMAX] && m_bounds[XMAX] >= box.m_bounds[XMIN] && m_bounds[YMIN] <= box.m_bounds[YMAX] && m_bounds[YMAX] >= box.m_bounds[YMIN] && m_bounds[ZMIN] <= box.m_bounds[ZMAX] && m_bounds[ZMAX] >= box.m_bounds[ZMIN];
        }
        
        inline bool intersects_ray( const real *origin, const real *direction ) const {
            assert( false );
            // TODO: implement this method
            return false;
        }
    };
    
};

#endif