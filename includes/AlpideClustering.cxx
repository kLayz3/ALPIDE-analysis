#include "AlpideClustering.h"
#include "AuxFunctions.h"
#include <cmath>

//#include<random>


using namespace std;
using namespace AlpideClustering;

uint32_t AlpideClustering::DistXY(const Point& p1, const Point& p2) { //obsolete
	return abs((int)p1.col - (int)p2.col) + abs((int)p1.row - (int)p2.row);
}

bool AlpideClustering::IsNeighbour(const Point& p1, const Point& p2) {
    if(abs((int)p1.col - (int)p2.col)<=1 && abs((int)p1.row - (int)p2.row)<=1) return true;
    return false;
}
bool AlpideClustering::IsInCluster(const Point& p0, const vector<Point>& cluster) {
    for(uint32_t i=0; i<cluster.size(); ++i) {
        if(IsNeighbour(p0, cluster[i]))
            return true;
    }
    return false;
}

vector<Point> AlpideClustering::MakeCluster(vector<Point>& hits) {
	vector<Point> cluster;
	Point p0 = hits[0];
	cluster.push_back(p0);

	/* If the 'distance' from a hits[i] to p0 is greater than 
	 * some value (set here as 100) then we can assume it's not part of the cluster.
	 * Hence it gets swapped with last 'interesting' hit, which 
	 * is at index j */

	int j = hits.size() - 1;

	bool pointFound = true;
	while(pointFound) {
		pointFound = false;
		for(int i=1; i<=j; ++i) {
			if(DistXY(p0, hits[i]) > 100) {
				QuickSwap(hits, i, j);
				--i; --j;
			}
			else if(IsInCluster(hits[i], cluster)) {
				cluster.push_back(hits[i]);
				QuickSwap(hits, i, j);
				QuickErase(hits, j);
				--i; --j;
				pointFound = true;
			}
		}
	}
	QuickErase(hits,0);
	return cluster;
}

vector<vector<Point>> AlpideClustering::ConstructClusters(uint32_t* ColArray, uint32_t* RowArray, uint32_t N, int veto) {
    vector<Point> hits;
    for(uint32_t i=0; i<N; ++i) 
        hits.push_back(Point(ColArray[i], RowArray[i]));

    vector<vector<Point>> clusters;

    /* Make isolated cluster around hits[0] with the MakeCluster(hits) call.
	 * Keep doing until hits size reaches 0
	 * Supply veto to not save clusters with size <= veto. Default 0 */

    while(hits.size() > 0) { /* Start a new cluster around hits[0] */
		auto cluster = MakeCluster(hits);
        if(cluster.size() > veto) clusters.push_back(cluster);
    }

    return clusters;
}

uint32_t AlpideClustering::FitCluster(const vector<Point>& cluster, double& uX, double& uY, double& sX, double& sY) {
    int N = cluster.size();
	if(N==0) {uX=0; uY=0; sX=0; sY=0; return 0;}
	uint32_t meanX{0}; uint64_t sigmaX{0}; //col
    uint32_t meanY{0}; uint64_t sigmaY{0}; //row
    for(auto p : cluster) {
        meanX += p.col; sigmaX += (uint64_t)(p.col*p.col);
        meanY += p.row; sigmaY += (uint64_t)(p.row*p.row);
    }
    uX  = (double)meanX/N;
    uY  = (double)meanY/N;
    
    /* For sigma=0 can underflow and throw exception, annoying */
    sX = (double)sigmaX/N - uX*uX;
    sX = (sX>0) ? sqrt(sX) : 0.0; 
    sY = (double)sigmaY/N - uY*uY;
    sY = (sY>0) ? sqrt(sY) : 0.0; 

    return N;
}

uint32_t AlpideClustering::FitCluster(const vector<Point>& cluster, double& uX, double& uY) {
    int N = cluster.size();
	if(N==0) {uX=0; uY=0; return 0;}
    uint32_t meanX(0); //col
    uint32_t meanY(0); //row
    for(auto& p : cluster) { 
        meanX += p.col; 
        meanY += p.row;
    }
    //Luke 
/*    std::random_device rd{};
    std::mt19937 gen{rd()};
   std::normal_distribution<> randX(0,26e-03/2.);
   std::normal_distribution<> randY(0,29e-03/2.);
  */
    uX  = ((double)meanX/N);//+(randX(gen)*(double)N);//Idea smeer wrt how many pixels x/2 as half the pixel width if we assume center of pixel? --May need modifying
    uY  = ((double)meanY/N);//+(randY(gen)*(double)N);
    return N;
}

