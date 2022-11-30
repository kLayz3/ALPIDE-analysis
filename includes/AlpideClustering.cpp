#include "AlpideClustering.h"
#include <cmath>
using namespace std;
using namespace AlpideClustering;

template<class T>  
void AlpideClustering::QuickErase(vector<T> &v, int i) {
    std::swap(v[i], v.back());
    v.pop_back();
}

bool AlpideClustering::IsNeighbour(const Point& p1, const Point& p2) {
    if(abs((int)p1.col - (int)p2.col)<=1 && abs((int)p1.row - (int)p2.row)<=1) return true;
    return false;
}
bool AlpideClustering::IsInCluster(const Point& p0, const vector<Point>& cluster) {
    for(unsigned i=0; i<cluster.size(); ++i) {
        if(IsNeighbour(p0, cluster[i]))
            return true;
    }   
    return false;
}

vector<Point> AlpideClustering::MakeCluster(Point& p0, vector<Point>& hits) {
    vector<Point> cluster;
    cluster.push_back(p0);
    bool pointFound(true);
    while(pointFound) {
        pointFound = false;
        for(unsigned i=0; i<hits.size(); ++i) {
            if(IsInCluster(hits[i], cluster)) {
                cluster.push_back(hits[i]);
                QuickErase(hits, i); 
                --i;
                pointFound = true;
            }   
        }   
    }   
    return cluster;
}

vector<vector<Point>> AlpideClustering::ConstructClusters(unsigned* ColArray, unsigned* RowArray, unsigned N, int veto) {
    vector<Point> hits;
    for(unsigned i=0; i<N; ++i) 
        hits.emplace_back(Point(ColArray[i], RowArray[i]));

    vector<vector<Point>> clusters;

    /* Initially mark a point an anchor. Go over all points left
     * in the initial vector. Find neighbours. Continue
     * the iteration until no new neighbours exist from
     * remaining points. That makes the cluster isolated. 
	 * Supply veto to not save clusters with size <= veto. Default 0 */

    while(hits.size() > 0) { /* Start a new cluster around hits[0] */
        Point p0 = hits[0];
        QuickErase(hits, 0); 
        auto cluster = MakeCluster(p0, hits);
        if(cluster.size() > veto) clusters.emplace_back(cluster); 
    }
    return clusters;
}

unsigned AlpideClustering::FitCluster(const vector<Point>& cluster, float& uX, float& uY, float& sX, float& sY) {
    int N = cluster.size();
	if(N==0) {uX=0; uY=0; sX=0; sY=0; return 0;}
    float meanX(0.); float sigmaX(0.); //col
    float meanY(0.); float sigmaY(0.); //row
    for(auto& p : cluster) {
        meanX += p.col; sigmaX += p.col*p.col;
        meanY += p.row; sigmaY += p.row*p.row;
    }
    meanX  = meanX/N;
    meanY  = meanY/N;
    
    /* For sigma=0 can underflow and throw exception, annoying */
    sigmaX = sigmaX/N - meanX*meanX;
    sigmaX = (sigmaX>0) ? sqrt(sigmaX) : 0; 
    sigmaY = sigmaY/N - meanY*meanY;
    sigmaY = (sigmaY>0) ? sqrt(sigmaY) : 0;

	/* Bind to arguments finally -> faster runtime when doing calculations on local stack variables 
	 * rather than pushing/popping passed-by-reference variables -- M.B. */
	uX = meanX; uY = meanY; sX = sigmaX; sY = sigmaY;
    return N;
}

unsigned AlpideClustering::FitCluster(const vector<Point>& cluster, float& uX, float& uY) {
    int N = cluster.size();
	if(N==0) {uX=0; uY=0; return 0;}
    float meanX(0.); //col
    float meanY(0.); //row
    for(auto& p : cluster) {
        meanX += p.col; 
        meanY += p.row;
    }
    meanX  = meanX/N;
    meanY  = meanY/N;	
	uX = meanX; uY = meanY;
    return N;
}
