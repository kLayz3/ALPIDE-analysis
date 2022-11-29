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
        Point(const std::pair<int,int> p) : col(p.first), row(p.second) {}; 
	    bool operator==(const Point& p) {
            if(col==p.col && row==p.row) return 1;
            return 0;
        }
		std::pair<int,int> toPair() {
			return std::make_pair(col, row);
		}
    };
    
    template<class T> void QuickErase(std::vector<T> &v, int i); 
    bool IsNeighbour(const Point& p1, const Point& p2);
    bool IsInCluster(const Point& p0, const std::vector<Point>& cluster);
    std::vector<Point> MakeCluster(Point& p0, std::vector<Point>& hits);
    std::vector<std::vector<Point>> ConstructClusters(unsigned* ColArray, unsigned* RowArray, unsigned N, int veto=0); 
    unsigned FitCluster(const std::vector<Point>& cluster, float& uX, float& uY, float& sX, float& sY);
}

#endif /* AlpideClustering.h */

