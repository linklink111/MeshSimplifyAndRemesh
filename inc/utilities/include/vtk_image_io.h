#ifndef VTK_IMAGE_IO_H
#define VTK_IMAGE_IO_H

#include<string>
#include<map>

#include<vtkSmartPointer.h>
#include<vtkFloatArray.h>
#include<vtkDoubleArray.h>
#include<vtkPointData.h>
#include<vtkCellData.h>
#include<vtkRectilinearGrid.h>
#include<vtkXMLRectilinearGridWriter.h>
#include<vtkXMLRectilinearGridReader.h>
#include<vtkImageData.h>
#include<vtkXMLImageDataWriter.h>
#include<vtkXMLImageDataReader.h>

namespace utilities {
		
	template< typename real, typename image >
	void load_vti( const char *filename, int *dim, real *aabb, std::map< std::string, image > &point_data, std::map< std::string, image > &cell_data ){
		vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
		reader->SetFileName( filename );
		reader->Update();
		
		vtkSmartPointer<vtkImageData> vti(reader->GetOutput());
		vti->GetDimensions(dim);
		dim[0]--; dim[1]--; dim[2]--;
		if( dim[2] <= 0 ) dim[2] = 1;
		
		int tdim[] = { dim[0]+1, dim[1]+1, dim[2]==1?dim[2]:dim[2]+1 };
		
		
		vtkSmartPointer<vtkPointData> pd(vti->GetPointData());
		for( int arr=0; arr<pd->GetNumberOfArrays(); arr++ ){
			vtkSmartPointer<vtkDoubleArray> tmp = vtkSmartPointer<vtkDoubleArray>::New();
			tmp->DeepCopy( pd->GetArray(arr) );
			tmp->SetName( pd->GetArray(arr)->GetName() );
			int nc = tmp->GetNumberOfComponents();
			image timage( tdim[0], tdim[1], tdim[2], nc );
			for( int i=0; i<tdim[0]; i++ ){
				for( int j=0; j<tdim[1]; j++ ){
					for( int k=0; k<tdim[2]; k++ ){
						int ijk[] = { i, j, k };
						for( int c=0; c<nc; c++ ){
							vtkIdType id = vti->ComputePointId( ijk );
							timage(i,j,k,c) = tmp->GetComponent( id, c );
						}
					}
				}
			}
			point_data[ std::string(tmp->GetName()) ] = timage;
		}
		vtkSmartPointer<vtkCellData> cd(vti->GetCellData());
		for( int arr=0; arr<cd->GetNumberOfArrays(); arr++ ){
			vtkSmartPointer<vtkDoubleArray> tmp = vtkSmartPointer<vtkDoubleArray>::New();
			tmp->DeepCopy( pd->GetArray(arr) );
			tmp->SetName( pd->GetArray(arr)->GetName() );
			int nc = tmp->GetNumberOfComponents();
			image timage( dim[0], dim[1], dim[2], nc );
			for( int i=0; i<dim[0]; i++ ){
				for( int j=0; j<dim[1]; j++ ){
					for( int k=0; k<dim[2]; k++ ){
						int ijk[] = { i, j, k };
						for( int c=0; c<nc; c++ ){
							vtkIdType id = vti->ComputeCellId( ijk );
							timage(i,j,k,c) = tmp->GetComponent( id, c );
						}
					}
				}
			}
			cell_data[ std::string(tmp->GetName()) ] = timage;
		}
	}
	
	template< typename real, typename image >
	void save_vti( const char *filename, const int *dim, const real *aabb, std::map< std::string, image > &point_data, std::map< std::string, image > &cell_data ){
		vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
		vtkSmartPointer<vtkImageData> vti = vtkSmartPointer<vtkImageData>::New();
		
		int tdim[] = { dim[0]+1, dim[1]+1, dim[2]==1?dim[2]:dim[2]+1 };
		vti->SetDimensions( tdim );
		
		for( typename std::map<std::string,image>::iterator iter=point_data.begin(); iter!=point_data.end(); iter++ ){
			image &img = iter->second;
			
			vtkSmartPointer<vtkFloatArray> vtkdata = vtkSmartPointer<vtkFloatArray>::New();
			vtkdata->SetName( iter->first.c_str() );
			vtkdata->SetNumberOfComponents( img.spectrum() );
			vtkdata->SetNumberOfTuples( tdim[0]*tdim[1]*tdim[2] );
			for( int i=0; i<tdim[0]; i++ ){
				for( int j=0; j<tdim[1]; j++ ){
					for( int k=0; k<tdim[2]; k++ ){
						int ijk[] = { i,j,k };
						vtkIdType id = vti->ComputePointId(ijk);
						for( int c=0; c<img.spectrum(); c++ ){
							vtkdata->SetComponent( id, c, img(i,j,k,c) );
						}
					}
				}
			}
			vti->GetPointData()->AddArray( vtkdata );
		}
		
		for( typename std::map<std::string,image>::iterator iter=cell_data.begin(); iter!=cell_data.end(); iter++ ){
			image &img = iter->second;
			
			vtkSmartPointer<vtkFloatArray> vtkdata = vtkSmartPointer<vtkFloatArray>::New();
			vtkdata->SetName( iter->first.c_str() );
			vtkdata->SetNumberOfComponents( img.spectrum() );
			vtkdata->SetNumberOfTuples( dim[0]*dim[1]*dim[2] );
			for( int i=0; i<dim[0]; i++ ){
				for( int j=0; j<dim[1]; j++ ){
					for( int k=0; k<dim[2]; k++ ){
						int ijk[] = { i,j,k };
						vtkIdType id = vti->ComputeCellId(ijk);
						for( int c=0; c<img.spectrum(); c++ ){
							vtkdata->SetComponent( id, c, img(i,j,k,c) );
						}
					}
				}
			}
			vti->GetCellData()->AddArray( vtkdata );
		}
		
		writer->SetInput( vti );
		writer->SetFileName( filename );
		writer->Write();
	}
	
	template< typename real, typename image >
	void load_vtr( const char *filename, int *dim, real *aabb, std::map< std::string, image > &point_data, std::map< std::string, image > &cell_data ){
		vtkSmartPointer<vtkXMLRectilinearGridReader> reader = vtkSmartPointer<vtkXMLRectilinearGridReader>::New();
		reader->SetFileName( filename );
		reader->Update();
		
		vtkSmartPointer<vtkRectilinearGrid> vtr(reader->GetOutput());
		vtr->GetDimensions(dim);
		dim[0]--; dim[1]--; dim[2]--;
		if( dim[2] <= 0 ) dim[2] = 1;
		
		int tdim[] = { dim[0]+1, dim[1]+1, dim[2]==1?dim[2]:dim[2]+1 };
		
		vtkSmartPointer<vtkDoubleArray> xcoords = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> ycoords = vtkSmartPointer<vtkDoubleArray>::New();
		vtkSmartPointer<vtkDoubleArray> zcoords = vtkSmartPointer<vtkDoubleArray>::New();
		xcoords->DeepCopy( vtr->GetXCoordinates() );
		ycoords->DeepCopy( vtr->GetYCoordinates() );
		zcoords->DeepCopy( vtr->GetZCoordinates() );
		aabb[0] = xcoords->GetValue(0); aabb[1] = xcoords->GetValue(dim[0]/*-1*/);
		aabb[2] = ycoords->GetValue(0); aabb[3] = ycoords->GetValue(dim[1]/*-1*/);
		aabb[4] = zcoords->GetValue(0); aabb[5] = zcoords->GetValue(dim[2]/*-1*/);
		
		vtkSmartPointer<vtkPointData> pd(vtr->GetPointData());
		for( int arr=0; arr<pd->GetNumberOfArrays(); arr++ ){
			vtkSmartPointer<vtkDoubleArray> tmp = vtkSmartPointer<vtkDoubleArray>::New();
			tmp->DeepCopy( pd->GetArray(arr) );
			tmp->SetName( pd->GetArray(arr)->GetName() );
			int nc = tmp->GetNumberOfComponents();
			image timage( tdim[0], tdim[1], tdim[2], nc );
			for( int i=0; i<tdim[0]; i++ ){
				for( int j=0; j<tdim[1]; j++ ){
					for( int k=0; k<tdim[2]; k++ ){
						int ijk[] = { i, j, k };
						for( int c=0; c<nc; c++ ){
							vtkIdType id = vtr->ComputePointId( ijk );
							timage(i,j,k,c) = tmp->GetComponent( id, c );
						}
					}
				}
			}
			point_data[ std::string(tmp->GetName()) ] = timage;
		}
		vtkSmartPointer<vtkCellData> cd(vtr->GetCellData());
		for( int arr=0; arr<cd->GetNumberOfArrays(); arr++ ){
			vtkSmartPointer<vtkDoubleArray> tmp = vtkSmartPointer<vtkDoubleArray>::New();
			tmp->DeepCopy( cd->GetArray(arr) );
			tmp->SetName( cd->GetArray(arr)->GetName() );
			int nc = tmp->GetNumberOfComponents();
			image timage( dim[0], dim[1], dim[2], nc );
			for( int i=0; i<dim[0]; i++ ){
				for( int j=0; j<dim[1]; j++ ){
					for( int k=0; k<dim[2]; k++ ){
						int ijk[] = { i, j, k };
						for( int c=0; c<nc; c++ ){
							vtkIdType id = vtr->ComputeCellId( ijk );
							timage(i,j,k,c) = tmp->GetComponent( id, c );
						}
					}
				}
			}
			cell_data[ std::string(tmp->GetName()) ] = timage;
		}
	}
	
	template< typename real, typename image >
	void save_vtr( const char *filename, const int *dim, const real *aabb, std::map< std::string, image > &point_data, std::map< std::string, image > &cell_data ){
		vtkSmartPointer<vtkXMLRectilinearGridWriter> writer = vtkSmartPointer<vtkXMLRectilinearGridWriter>::New();
		vtkSmartPointer<vtkRectilinearGrid> vtr = vtkSmartPointer<vtkRectilinearGrid>::New();
		vtkSmartPointer<vtkFloatArray> xcoords = vtkSmartPointer<vtkFloatArray>::New();
		vtkSmartPointer<vtkFloatArray> ycoords = vtkSmartPointer<vtkFloatArray>::New();
		vtkSmartPointer<vtkFloatArray> zcoords = vtkSmartPointer<vtkFloatArray>::New();
		real h[] = { (aabb[1]-aabb[0])/real(dim[0]), (aabb[3]-aabb[2])/real(dim[1]), (aabb[5]-aabb[4])/real(dim[2]) };
		for( int i=0; i<dim[0]+1; i++ ){
			xcoords->InsertNextValue( aabb[0]+i*h[0] );
		}
		for( int i=0; i<dim[1]+1; i++ ){
			ycoords->InsertNextValue( aabb[2]+i*h[1] );
		}
		if( dim[2] > 1 ){
			for( int i=0; i<dim[2]+1; i++ ){
				zcoords->InsertNextValue( aabb[4]+i*h[2] );
			}
		} else {
			zcoords->InsertNextValue( 0.0 );
		}
		int tdim[] = { dim[0]+1, dim[1]+1, dim[2]==1?dim[2]:dim[2]+1 };
		vtr->SetDimensions( tdim );
		vtr->SetXCoordinates( xcoords );
		vtr->SetYCoordinates( ycoords );
		vtr->SetZCoordinates( zcoords );
		
		for( typename std::map<std::string,image>::iterator iter=point_data.begin(); iter!=point_data.end(); iter++ ){
			image &img = iter->second;
			
			vtkSmartPointer<vtkFloatArray> vtkdata = vtkSmartPointer<vtkFloatArray>::New();
			vtkdata->SetName( iter->first.c_str() );
			vtkdata->SetNumberOfComponents( img.spectrum() );
			vtkdata->SetNumberOfTuples( tdim[0]*tdim[1]*tdim[2] );
			for( int i=0; i<tdim[0]; i++ ){
				for( int j=0; j<tdim[1]; j++ ){
					for( int k=0; k<tdim[2]; k++ ){
						int ijk[] = { i,j,k };
						vtkIdType id = vtr->ComputePointId(ijk);
						for( int c=0; c<img.spectrum(); c++ ){
							vtkdata->SetComponent( id, c, img(i,j,k,c) );
						}
					}
				}
			}
			vtr->GetPointData()->AddArray( vtkdata );
		}
		
		for( typename std::map<std::string,image>::iterator iter=cell_data.begin(); iter!=cell_data.end(); iter++ ){
			image &img = iter->second;
			
			vtkSmartPointer<vtkFloatArray> vtkdata = vtkSmartPointer<vtkFloatArray>::New();
			vtkdata->SetName( iter->first.c_str() );
			vtkdata->SetNumberOfComponents( img.spectrum() );
			vtkdata->SetNumberOfTuples( dim[0]*dim[1]*dim[2] );
			for( int i=0; i<dim[0]; i++ ){
				for( int j=0; j<dim[1]; j++ ){
					for( int k=0; k<dim[2]; k++ ){
						int ijk[] = { i,j,k };
						vtkIdType id = vtr->ComputeCellId(ijk);
						for( int c=0; c<img.spectrum(); c++ ){
							vtkdata->SetComponent( id, c, img(i,j,k,c) );
						}
					}
				}
			}
			vtr->GetCellData()->AddArray( vtkdata );
		}
		
		writer->SetInput( vtr );
		writer->SetFileName( filename );
		writer->Write();
	}
	
};

#endif
