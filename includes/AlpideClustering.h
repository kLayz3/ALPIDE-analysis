#ifndef ALPIDECLUSTERING_H
#define ALPIDECLUSTERING_H
#include <vector>
#include <tuple>

namespace AlpideClustering {
    struct Point {
        unsigned col;
        unsigned row;
        Point() = default;
        Point(unsigned col, unsigned row): col(col), row(row) {}; 
        Point(const Point& p) : col(p.col), row(p.row) {}; 
        bool operator==(const Point& p) {
            if(col==p1.col && row==p1.row) return 1;
            return 0;
        }
    };
    
    template<class T> void QuickErase(std::vector<T> &v, int i); 
    bool IsNeighbour(const Point& p1, const Point& p2);
    bool IsInCluster(const Point& p0, const std::vector<Point>& cluster);
    std::vector<Point> MakeCluster(Point& p0, std::vector<Point>& hits);
    std::vector<std::vector<Point>> ConstructClusters(unsigned* ColArray, unsigned* RowArray, unsigned N); 
    std::tuple<double,double,double,double,int> FitCluster(const std::vector<Point>& cluster);
}

#endif /* AlpideClustering.h */

