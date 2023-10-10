#ifndef ALPIDE_CLUSTERING_H
#define ALPIDE_CLUSTERING_H

#include <vector>
#include <tuple>
#include <cstdint>

namespace AlpideClustering {
    struct Point {
        uint32_t col;
        uint32_t row;
        Point() = default;
        Point(uint32_t col, uint32_t row): col(col), row(row) {}; 
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

	typedef std::vector<Point> Cluster;
	typedef std::vector<Cluster> ClusterVec;

	uint32_t DistXY(const Point p1, const Point p2);
    inline bool IsNeighbour(const Point p1, const Point p2);
    bool IsInCluster(const Point p0, const std::vector<Point>& cluster);
    std::vector<Point> MakeCluster(std::vector<Point>& hits);

    std::vector<std::pair<uint32_t, ClusterVec>> ConstructClusters(uint32_t* chipArray, uint32_t* ColArray, uint32_t* RowArray, uint32_t N, int veto=0); 
	// ----> chipId ------^^^^^^^^
    uint32_t FitCluster(const std::vector<Point>& cluster, double& uX, double& uY, double& sX, double& sY);
    uint32_t FitCluster(const std::vector<Point>& cluster, double& uX, double& uY);

}

#endif 
