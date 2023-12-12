#pragma once
#include <vtkRenderWindow.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <vtkScalarBarRepresentation.h>
#include <vtkScalarBarWidget.h>
#include <vtkScalarBarActor.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkCubeAxesActor.h>
struct Points_xyz_range
{
	Points_xyz_range() :x_min(NAN), x_max(NAN),y_min(NAN), y_max(NAN),z_min(NAN), z_max(NAN)
	{
	}
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double z_min;
	double z_max;
};
class VTK_plot
{
public:
	static void set_lut_mode(vtkSmartPointer<vtkLookupTable> vtk_lut, int color_mode);
	static void init_lut(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, vtkSmartPointer<vtkLookupTable> vtk_lut, vtkSmartPointer<vtkScalarBarActor> scalarBarActor_for_vtk, vtkSmartPointer<vtkScalarBarWidget> scalarbar_Widget_member_for_vtk, vtkSmartPointer<vtkScalarBarRepresentation> scalarBarRep_for_vtk);
	static Points_xyz_range get_xyz_range(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud);
	static Points_xyz_range get_vec_range(std::vector<double> vec);
	static void set_clouds_color_other(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, std::vector<double> other_vec,
		vtkSmartPointer<vtkLookupTable> vtk_lut, Points_xyz_range o_range);
	static void set_clouds_color_xyz(pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud, vtkSmartPointer<vtkLookupTable> vtk_lut, Points_xyz_range xyz_range, int color_mode);
	static void set_grid_mode(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, vtkSmartPointer<vtkCubeAxesActor> cubeAxesActor, Points_xyz_range xyz_range, bool show_mode);
	static void set_lut_show(boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer, vtkSmartPointer<vtkLookupTable> vtk_lut, vtkSmartPointer<vtkScalarBarActor> scalarBarActor_for_vtk, double min_val,double max_val, std::string label);
};

