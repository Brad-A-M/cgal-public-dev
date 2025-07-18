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



const double minDistance = .01;
 
double weight (Point p, std::vector<Point> neighbors)
{
  double weight = 0;
  for(Point n : neighbors)
  {
      
      weight = weight + pow(1-(sqrt(CGAL::squared_distance(p,n))/minDistance),8);
  }
    
  return weight;
}

class BidirectionalMap {
public:
    // Maps to store the point-to-int and int-to-point mappings
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

void update_indices(std::vector<Point> neighbors, std::vector<Weighted> data, BidirectionalMap bimap){
    //traverse neighbors
    //  - check their children in data, 2i + 1 and 2i + 2, move them down heap while one child is bigger 
    //  - reassign the index in biMap 

    // std::cout << "Before traversing neighbors within update_indices\n" << std::endl;

    for(Point n : neighbors){
        double curr_w = 0;
        double child1_w = 0;
        double child2_w = 0;

         std::cout << "Right before while loop within update_indices\n" << std::endl;

        do{
            int n_index = bimap.getInt(n);

            if(n_index >= (data.size() / 2) - 1){
                break;
            }

            Weighted curr = data.at(n_index);
            Weighted child1 = data.at(2 * n_index + 1);
            Weighted child2 = data.at(2 * n_index + 2);

            curr_w = curr.weight;
            child1_w = child1.weight;
            child2_w = child2.weight;
            std::cout<< "current weight: "<< curr.weight << " c1 weight: " << child1.weight << " c2 weight: " << child2.weight << std::endl;

            if(child1.weight > child2.weight){
                if(child1.weight > curr.weight){
                    data[n_index] = child1;
                    data[2 * n_index + 1] = curr;

                    // std::cout << "Beofre swap w/ child1\n" << std::endl;

                    bimap.swap(n, child1.point);
                }
            }else if(child2.weight > child1.weight){
                if(child2.weight > curr.weight){
                    data[n_index] = child2;
                    data[2 * n_index + 2] = curr;

                    // std::cout << "Beofre swap w/ child2\n" << std::endl;

                    bimap.swap(n, child2.point);
                }
            }
        }while(child1_w > curr_w || child2_w > curr_w);

        // std::cout << "After while loop\n" << std::endl;
    }
    // std::cout << "Done in update indices\n";
}
 
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
                              CGAL::parameters::number_of_points_per_area_unit(20000));
   
   
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
  Fuzzy_circle default_range(query, minDistance);
    
    
  std::vector<Point> result;
  tree.search(std::back_inserter(result), default_range);
    
  std::cout << "\nPoints in circle with center: " << query << " and radius:"<< minDistance << std::endl;
    
 
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
      Fuzzy_circle default_range(p,minDistance);
      tree.search(std::back_inserter(neighbors), default_range);
      double w = weight(p,neighbors);
     // weightedHeap.push({p,w,neighbors});
      data.push_back({p,w,neighbors});
      neighbors.clear();
      
  }
    
   
  std::cout << "First weight before heap: "<< data.begin()->weight << "?" << std::endl;
  std::make_heap(data.begin(), data.end(),CompareWeighted());
  std::cout << "First weight after heap: "<< data.begin()->weight << "?" << std::endl;
   
    /*
  std::vector<Point> max_weight;
  max_weight.push_back(data.begin()->point);
  std::ofstream out3("max_weight.xyz");
  out3 << std::setprecision(17);
  std::copy(max_weight.begin(), max_weight.end(), std::ostream_iterator<Point>(out3, "\n"));
  out3.close();
*/
    
 BidirectionalMap bimap;
 for(int i = 0; i < data.size(); ++i)
 {
     bimap.insert(data[i].point,i);
 }
    /*Commented shows an example of bimap
        std::cout << "Point -> " << data[42].point << " has index -> " << bimap.getInt(data[42].point) << std::endl;
    
      std::cout << "Index -> " << 42 << " has point -> " << bimap.getPnt(42) << std::endl;
*/
    /*
 for(int i = 0; i < 10000; i++){
    //traverse data vector from 0 -> k
    //for heaviest weight 
    // - traverse through neighbors and update weights
    for(Point p : data[i].neighbors){
        int ind = bimap.getInt(p);
        data[ind].weight = data[ind].weight - pow((1-(sqrt(CGAL::squared_distance(p,data[i].point))/minDistance)),8);
    }

    update_indices(data[i].neighbors, data, bimap);

    std::cout << i << std::endl;
 }

   
 std::vector<Point> sample;
 //std::cout << "Points in the sample: \n" << std::endl;
 for(int i = 10000; i < data.size(); i++){
    //collect points and print them out
    sample.push_back(data[i].point);
   // std::cout << data[i].point << "\n" << std::endl;
 }
    
    
 std::ofstream out2("finial_sample.xyz");
 out2 << std::setprecision(17);
 std::copy(sample.begin(), sample.end(), std::ostream_iterator<Point>(out2, "\n"));
 out2.close();
*/
    
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
