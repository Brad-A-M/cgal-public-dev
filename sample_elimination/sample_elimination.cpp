#include <CGAL/Simple_cartesian.h>
#include <CGAL/point_generators_3.h>
//#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Search_traits_3.h>
//#include <CGAL/Orthogonal_incremental_neighbor_search.h>

 
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
 
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Point_set_3.h>
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <queue>
#include <unordered_map>

#include <CGAL/Kd_tree.h>
#include <CGAL/algorithm.h>
#include <CGAL/Fuzzy_sphere.h>
 
typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef K::Point_3                                            Point;
 
typedef CGAL::Point_set_3<Point>                              Point_set;

 
typedef CGAL::Surface_mesh<Point>                             Mesh;
typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;
typedef boost::graph_traits<Mesh>::face_descriptor            face_descriptor;

typedef CGAL::Search_traits_3<K> Traits;
typedef CGAL::Fuzzy_sphere<Traits> Fuzzy_circle;
typedef CGAL::Kd_tree<Traits> Tree;


namespace PMP = CGAL::Polygon_mesh_processing;

struct Weighted {
    Point point;
    double weight;
    std::vector<Point> neighbors;

    Weighted(Point point, double weight, std::vector<Point> neighbors) : point(point), weight(weight), neighbors(neighbors) {}
   
};

struct CompareWeighted {
    bool operator()(const Weighted& a, const Weighted& b) {
        return a.weight < b.weight; // Max-heap
    }
};



const double minDistance = .02;
 
double weight (Point p, std::vector<Point> neighbors)
{
  double weight = 0;
  for(Point n : neighbors)
  {
      weight = weight + pow((1-sqrt(CGAL::squared_distance(p,n)))/(minDistance),8);
  }
    
  return weight;
}

void update_indices(std::vector<Point> neighbors, std::vector<Weighted> data, BidirectionalMap bimap){
    //traverse neighbors
    //  - check their children in data, 2i + 1 and 2i + 2, move them down heap while one child is bigger 
    //  - reassign the index in biMap 

    for(Point n : neighbors){
        do{
            n_index = bimap.getInt(n);

            if(n_index >= (data.size() / 2) - 1){
                break;
            }
            
            Point child1 = data.at(2 * n_index + 1);
            Point child2 = data.at(2 * n_index + 2);

            if(child1.weight > child2.weight){
                if(child1.weight > n.weight){
                    data[n_index] = child1;
                    data[2 * n_index + 1] = n;

                    bimap.swap(n, child1);
                }
            }else if(child2.weight > child1.weight){
                if(child2.weight > n.weight){
                    data[n_index] = child2;
                    data[2 * n_index + 2] = n;

                    bimap.swap(n, child2);
                }
            }
        }while(child1.weight > n.weight || child2.weight > n)
    }
}

class BidirectionalMap {
public:
    // Maps to store the string-to-int and int-to-string mappings
    std::unordered_map<Point, int> pntToInt;
    std::unordered_map<int, Point> intToPnt;

    // Insert a pair into the map
    void insert(const Point& pnt, int num) {
        pntToInt[pnt] = num;
        intToPnt[num] = pnt;
    }

    //Swap two points indices in the map
    void swap(const Point& pnt1, const Point& pnt2){
        int pnt1_ind = getInt(pnt1);
        int pnt2_ind = getInt(pnt2);

        insert(pnt1, pnt2_ind);
        insert(pnt2, pnt1_ind);
    }

    // Get the integer corresponding to the point
    int getInt(const Point& pnt) {
        if (pntToInt.find(pnt) != pntToInt.end()) {
            return pntToInt[pnt];
        }
        throw std::invalid_argument("Point not found!");
    }

    // Get the point corresponding to the integer
    Point getPnt(int num) {
        if (intToPnt.find(num) != intToPnt.end()) {
            return intToPnt[num];
        }
        throw std::invalid_argument("Integer not found!");
    }
};


 
int main(int argc, char* argv[])
{
    
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/eight.off");
   
  Mesh mesh;
  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
     std::cerr << "Invalid input." << std::endl;
     return 1;
  }
   
  const int points_per_face = (argc > 2) ? std::stoi(argv[2]) : 30;
   
  std::vector<Point> points;
  PMP::sample_triangle_mesh(mesh,
                              std::back_inserter(points),
                              CGAL::parameters::number_of_points_per_area_unit(50000));
   
   
  //Point_set point_set;
  //PMP::sample_triangle_mesh(mesh,
                           //   point_set.point_back_inserter());
   
  std::cout << "Initial sample size: " << points.size() << std::endl;
    
  std::cout << "Number of faces in my mesh: " << mesh.number_of_faces() << std::endl;
  std::ofstream out("initial_sample.xyz");
  out << std::setprecision(17);
  std::copy(points.begin(), points.end(), std::ostream_iterator<Point>(out, "\n"));
  out.close();
    
    

 
  // Build tree in parallel
  Tree tree(points.begin(), points.end());
  tree.build<CGAL::Parallel_tag>();
    
  Point query = points.front();
  Fuzzy_circle default_range(query, .02);
    
    
  std::vector<Point> result;
  tree.search(std::back_inserter(result), default_range);
    
  std::cout << "\nPoints in circle with center: " << query << " and radius: 0.02" << std::endl;
    
 
  for (size_t i = 0; i < result.size(); ++i) {
            std::cout << result[i] << "\n ";
        }
  
    
  std::ofstream out1("ball_o_points.xyz");
  out1 << std::setprecision(17);
  std::copy(result.begin(), result.end(), std::ostream_iterator<Point>(out1, "\n"));
  out1.close();
    
  //Building my data structures
    
    
  std::vector<Point> neighbors;

  std::vector<Weighted> data;
   
  for (Point p : points)
  {
      Fuzzy_circle default_range(p,.02);
      tree.search(std::back_inserter(neighbors), default_range);
      double w = weight(p,neighbors);
     // weightedHeap.push({p,w,neighbors});
      data.push_back({p,w,neighbors});
      neighbors.clear();
      
  }
 
  std::cout << "First weight before heap: "<< data.begin()->weight << "?" << std::endl;
  std::make_heap(data.begin(), data.end(),CompareWeighted());
  std::cout << "First weight after heap: "<< data.begin()->weight << "?" << std::endl;
    

 BidirectionalMap bimap;
 for(int i = 0; i < data.size(); ++i)
 {
     bimap.insert(data[i].point,i);
 }
  
 std::cout << "Point -> " << data[42].point << " has point -> " << bimap.getInt(data[42].point) << std::endl;
    
 std::cout << "Index -> " << 42 << " has point -> " << bimap.getPnt(42) << std::endl;


 for(int i = 0; i < 40000; i++){
    //traverse data vector from 0 -> k
    //for heaviest weight 
    // - traverse through neighbors and update weights

    //make function for both of the below:
    // - move adjusted neighbors down heap
    // - update indices
 }

 for(int i = 40000; i < data.size(), i++){
    //collect points and print them out
 }
    
/*
 
  // Query tree in parallel
  std::vector<std::vector<Point> > neighbors (points.size());
  tbb::parallel_for (tbb::blocked_range<std::size_t> (0, points.size()),
                     [&](const tbb::blocked_range<std::size_t>& r)
                     {
                       for (std::size_t s = r.begin(); s != r.end(); ++ s)
                       {
                         // Neighbor search can be instantiated from
                         // several threads at the same time
                         Neighbor_search search (tree, points[s], k);
                         neighbors[s].reserve(k);
 
                         // neighbor search returns a set of pair of
                         // point and distance <Point_3,FT>, here we
                         // keep the points only
                         for (const Point_with_distance pwd : search)
                           neighbors[s].push_back (pwd.first);
                       }
                     });
 */
  return 0;
}
