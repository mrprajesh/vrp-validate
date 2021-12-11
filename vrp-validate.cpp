 /*
 0* main.cpp
 *
 * Copyright 2020 Team Resilience of IIT Madras for 12th DIMACS Contest
 *
 * Rajesh Pandian M <mrpaj at gmail dot com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 *
 *
 *
 * vr-validate.cpp
 */

// Changelog
// BEGIN code                             Fri, 28-May-2021, 12:42:12 IST
// TODO VALIDATE                          Fri, 28-May-2021, 12:42:36 IST
// DONE VALIDATE                          ...
// ADDD BUGFix on float inputs            Thu, 09-Dec-2021, 12:12:07 IST


#include <iostream>   // also has pair!?
#include <vector>
#include <sstream>    //istringstream
#include <algorithm>  //transform
#include <cmath>      //sqrt
#include <map>
#include <set>
#include <list>
#include <climits>
#include <cstring>

#include <chrono> //timing CPU

#define DEBUG if(DEBUGCODE)
#define DEFAULT -23432

unsigned DEBUGCODE = 0;

// Let's have ~~int~~ Ok Templated
template<typename T>
class CVRPInOut{
public:
  CVRPInOut(){
    //! routes.reserve(n);
  }
  //! CVRPInOut(){}
  ~CVRPInOut(){}
  void addPoint(T x, T y){
    pts.push_back({x,y});
  }
  void addDemand(T d){
    demand.push_back(d);
  }

  void setCapacity(T cap){
    residueCap = cap;
    vCapacity = cap;
  }
  std::pair<T,T> getPoint(T id){
    return pts[id-1];
  }

  bool addRouteVertex(T v){
    if(residueCap - demand[v] >= 0){
      aRoute.push_back(v);
      //! routes[curRoute].push_back(v);
      residueCap = residueCap - demand[v];
    }
    else{   //
      routes.push_back(aRoute);
      aRoute.clear();
      aRoute.push_back(v);
      //! routes[++curRoute].push_back(v);
      residueCap = vCapacity - demand[v];
    }


    return true;
  }
  std::vector<T> aRoute;

  std::vector< std::list<T>> routes;
  std::vector< std::vector<T>> routesVector;

  std::vector<T> allNodes;

  float cost;

  std::vector<T>  rDemand;  //Route's commulative demand
  std::vector<T>  r;        //Route's id

  T residueCap;

  T n;
  T vCapacity;
  std::string     dim;
  std::vector<T>  demand;

  std::vector<std::pair<T,T>> pts;
  private:
  T curRoute = 0;
  //! later
  //! operator () {
  //! }
};

template<typename T>
class cvrpOut{
public:
  cvrpOut(){
    vCapacity  = -1;
    residueCap = -1;
  }

  ~cvrpOut(){}

  cvrpOut(T c){
    vCapacity = c;
    residueCap = c;
  }
  //! bool addRouteVertex(T r){
    //! if(residueCap  )
    //! aRoute.push_back(r);
    //! return true;
  //! }
  std::vector<T> aRoute;
  std::vector< std::vector<T>> routes;
  std::vector< std::list<T>> r;
  std::vector<T> rDemand;

  T residueCap;
  T vCapacity;
};

template<typename T>
void extract(std::string &line, T &val){
  std::string dummy; // for code and ":"
  std::istringstream iss(line);
  iss >> dummy >> dummy >> val;
  //! std::cout<< dummy << " val " << val << '\n';
}

template<typename T>
void extractPoints(CVRPInOut<T> &inputs){
  T id, x, y;
  std::string line, xStr, yStr;
  for(T ii=1, end=inputs.n; ii <= end && std::getline(std::cin, line) ; ++ii){
    std::istringstream iss(line);

    iss >> id >> xStr >> yStr;    // BUG FIX
    x=stof(xStr);
    y=stof(yStr);
    //! std::cout<< id  << " " << x << ","<< y <<" " << '\n';

    if(id != ii) {
      std::cout<< "WARN: Format maybe incorrect!" << '\n';
      return;
    }
    inputs.addPoint(x,y); //inserting pts 0 to n-1  and NO check of id while inserting. So, if block will WARN
  }
}
template<typename T>
void extractDemands(CVRPInOut<T> &inputs){
  T id, d;
  std::string line;
  for(T ii=1, end=inputs.n; ii <= end && std::getline(std::cin, line) ; ++ii){
    std::istringstream iss(line);
    iss >> id >> d;
    if(id != ii) {
      std::cout<< "WARN: Format maybe incorrect!" << '\n';
      return;
    }
    inputs.addDemand(d);
  }
}


template<typename T>
void readVRP(CVRPInOut<T> &inputs){
  std::string line, code, dummy;
  T depot = DEFAULT;
  T val = DEFAULT;

  while (std::getline(std::cin, line)) {
    std::istringstream iss(line);
    //! std::cout<< line << '\n';
    if(iss >> code && code != ""){
      //! std::cout<< code << '\n';
      transform(code.begin(), code.end(), code.begin(), ::toupper); // just in case if not upper
      if(code == "DIMENSION"){
        extract<T>(line, inputs.n);
        //! std::cout<< "n:" << inputs.n  << '\n';
      }else if (code == "CAPACITY") {
        extract<T>(line, val);
        inputs.setCapacity(val);
      }else if (code == "EDGE_WEIGHT_TYPE"){
        extract<std::string>(line, inputs.dim); //store 2 or 2D; let's have as str
      }else if (code == "NODE_COORD_SECTION"){
        extractPoints<T>(inputs);
      }else if (code == "DEMAND_SECTION"){
        extractDemands<T>(inputs);
      }else if (code == "DEPOT_SECTION"){
        std::cin >> depot >> dummy;
      }else if (code == "EOF"){
        if(depot != 1 || depot == DEFAULT){
          std::cout<< "WARN: Depot may not be 1" << '\n';
          exit(0);
        }
        DEBUG std::cout<< "INPUT - Done" << '\n';
        return;
      }
      else if (code == "NAME" || code =="COMMENT" || code == "TYPE"  ){
        //do nothing
      }
      else{
        std::cout << "WARN: CODE "<< code << '\n';
      }
    }
  }
}

float euDist(int x1, int y1, int x2, int y2){
  return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

template <typename T1, typename T2>
void calDist(const CVRPInOut<T1> &inputs, std::vector< std::vector<T2>> &dist){
  T1 size=inputs.n;
  dist.resize(size);
  for(T1 ii=0 ; ii < size; ++ii)
    dist[ii].resize(size,0);

  for(T1 ii=0 ; ii < size; ++ii){  // ROOM for parallelism
    dist[ii].resize(size,0);
    for(T1 jj=0 ; jj < size; ++jj){
      if(ii<jj){  //symmetric
        T1 x1 = inputs.pts[ii].first ;
        T1 y1 = inputs.pts[ii].second;
        T1 x2 = inputs.pts[jj].first ;
        T1 y2 = inputs.pts[jj].second;
        T2 distance = euDist(x1,y1,x2,y2);
        DEBUG std::cout<< distance << '\n';
        dist[ii][jj] = distance;
        dist[jj][ii] = distance;
      }
    }
  }
}


template <typename T1, typename T2>
void computeSavings(const CVRPInOut<T1> &inputs, std::vector< std::vector<T2>> &dist,
  //! std::set<std::pair<T2, std::pair<T1,T1>> , std::greater<std::pair<T2, std::pair<T1,T1>>>  > &savings
  auto &savings
){
  T1 size=inputs.n;
  const T1 DEPOT = 0;
  for(T1 ii=0 ; ii < size; ++ii){  // ROOM for parallelism
    dist[ii].resize(size,0);
    for(T1 jj=0 ; jj < size; ++jj){
      if(ii<jj){  //symmetric
        float sij= dist[DEPOT][ii] + dist[DEPOT][jj] - dist[ii][jj];
        savings.insert({sij,{ii,jj}});
      }
    }
  }

  DEBUG for(auto a: savings)
    std::cout<< a.first << " " << a.second.first << " " << a.second.second  << '\n';
}

template <typename T1>
void printRoutes(const CVRPInOut<T1> &io){
  T1 i = 0;
  std::cout<< "==========\n";
  for(auto &aRoute: io.routes){
    std::cout<< i++ <<": " << ' ';
    for(auto &node: aRoute)
      std::cout<< node <<' ';
    std::cout<< '\n';
  }
  std::cout<< "---\n";
  i = 0;
  for (auto &cost : io.rDemand)
    std::cout<< i++ << " "<< cost << '\n';
  std::cout<< "---\n";
  i = 0;
  for (auto &a : io.r)
    std::cout<< " r["<< i++ << "] = "<< a;
  std::cout<< "" << '\n';
}

template <typename T1>
bool addRoute(CVRPInOut<T1> &io, T1 i, T1 j){
  std::list <T1> aRoute;
  aRoute.push_back(i);
  aRoute.push_back(j);
  T1 rCost = io.demand[i]+io.demand[j];
  if( rCost <= io.vCapacity){
    io.routes.push_back(aRoute);  // new routes created
    io.rDemand.push_back(rCost);
    io.r[i] = io.r[j] = io.routes.size() - 1;
    return true;
  }
  DEBUG std::cout<< "EXCEEDED Capacity " ;
  return false;
}

template <typename T1>
bool addToRoute(CVRPInOut<T1> &io, T1 i, T1 j){
  T1 ri = io.r[i];
  //! T1 rj = io.r[j];
  T1 rCost = io.rDemand[ri] + io.demand[j]; //mistake?
  if (rCost <= io.vCapacity ){
    //! T1 f = io.routes[ri].front();
    //! T1 b = io.routes[ri].back() ;

    //! std::cout<< "F: "<< f << " B:"<< b << " Ri: " << ri << " Rj:" << rj << ' ';
    if( io.routes[ri].front() == i){
      io.r[j] = ri;
      //! std::cout<< "@  j: " << j << " R[j]:" << io.r[j] << ' ';
      io.rDemand[ri] += io.demand[j];
      io.routes[ri].push_front(j);
      return true;
    }else if (io.routes[ri].back() == i) {
      io.r[j] = ri;
      io.rDemand[ri] += io.demand[j];
      io.routes[ri].push_back(j);
      return true;
    }
    else {
      DEBUG std::cout<< "NOT INTERIOR " ;
      return false;
    }
  }
  DEBUG std::cout<< "EXCEEDED Capacity " ;
  return false;
}

template <typename T1>
bool mergeRoute(CVRPInOut<T1> &io, T1 i, T1 j){
  T1 ri = io.r[i];
  T1 rj = io.r[j];

  if (ri == rj){ // ==
    DEBUG std::cout<< "NOT Suited - SAME BAG1 " ;
    return false;
  }
  T1 rCost = io.rDemand[ri] + io.rDemand[rj];
  if(rCost <= io.vCapacity){
    if(ri < rj) {
      DEBUG std::cout<< " 1 " << " - ";
      io.routes[ri].merge(io.routes[rj]);
      io.rDemand[ri] +=io.rDemand[rj];
      for(auto &v: io.routes[rj]){
        io.r[v] = ri;
      }
      io.routes[rj].clear();
      io.rDemand[rj] = 0;
      return true;
    }
    else if (ri > rj){
      DEBUG std::cout<< " 2 " << " - ";
      io.routes[rj].merge(io.routes[ri]);
      io.rDemand[rj] +=io.rDemand[ri];
      for(auto &v: io.routes[ri]){
        io.r[v] = rj;
      }
      io.routes[ri].clear();
      io.rDemand[ri] = 0;
      return true;
    }
    else{ // ==
      DEBUG std::cout<< "NOT Suited - SAME BAG2 " ;
      return false;
    }
  }
  DEBUG std::cout<< "EXCEEDED Capacity " ;
  return false;
}

template <typename T1, typename T2>
void procesSavingsList(CVRPInOut<T1> &io, const std::vector< std::vector<T2>> &dist, auto &savings
){
  auto n = io.n;
  const auto NOTROUTED = -1;
  //! -1: NOTROUTED //
  //!  0: DEPOT     // Depot is routed in all the routes. So a special tag
  //! >1: ROUTED    // Routes starts from 1

  io.r.resize(n, NOTROUTED);

  // INITS TO MAKE INDEX BEGIN FROM 1
  std::list<T1> aList; aList.push_back(0);
  io.routes.push_back(aList);
  io.rDemand.push_back(0);

  DEBUG std::cout<< "rSize "<< io.r.size() << '\n';

  io.r[0] = 0;

  // savings = funct()

  for (auto &a : savings){

    T1 i = a.second.first;
    T1 j = a.second.second;

    bool isRoutedRi = (io.r[i] > 0) ? true: false;
    bool isRoutedRj = (io.r[j] > 0) ? true: false;

    if(i == 0 || j == 0){
      DEBUG std::cout<< "WARN depoT" << '\n'; // TESTINGg
      break;
    }

    if (!isRoutedRi && !isRoutedRj){
      DEBUG std::cout<< "CASE 1 (" <<  i <<","<< j << ") ";

      if(addRoute<T1>(io, i, j)){
        DEBUG std::cout<< "Success" << '\n';
      }else {
        DEBUG std::cout<< "Failure" << '\n';
      }
    }
    else if (isRoutedRi && !isRoutedRj){
      DEBUG std::cout<< "CASE 2.1L (*" <<  i <<","<< j << ") ";
      if(addToRoute(io, i,j)) {
        DEBUG std::cout<< "Success" << '\n';
      }else {
        DEBUG std::cout<< "Failure" << '\n';
      }
    }
    else if(!isRoutedRi && isRoutedRj) {
      DEBUG std::cout<< "CASE 2.2R (" <<  i <<","<< j << "*) ";
      if(addToRoute(io, j,i)) {
        DEBUG std::cout<< "Success" << '\n';
      }else{
        DEBUG std::cout<< "Failure" << '\n';
      }

    }
    else if (isRoutedRi && isRoutedRj){
      DEBUG std::cout<< "CASE 3 (" <<  i <<","<< j << ") ";//Different route and not interior
      if(mergeRoute(io, i,j)) {
        DEBUG std::cout<< "Success" << '\n';
      }else{
        DEBUG std::cout<< "Failure" << '\n';
      }

    }
    DEBUG printRoutes<T1>(io);
  }
  DEBUG std::cout<< "Total Routes: "<< io.routes.size() << '\n';
}

template<typename T1,typename T2>
void printDist(const std::vector< std::vector<T2>> &dist){
  T1 size = dist.size();
  for(T1 ii=0 ; ii < size; ++ii){
    for(T1 jj=0 ; jj < size; ++jj){
      DEBUG std::cout<< dist[ii][jj] << "\t" ;
    }
    DEBUG std::cout<< " " << '\n';
  }
}

template<typename T>
void printInputs(const CVRPInOut<T> &inputs){
  std::cout<< "N\t:"    << inputs.n           << '\n';
  std::cout<< "DIM\t:"  << inputs.dim         << '\n';
  std::cout<< "CAP\t:"  << inputs.vCapacity   << '\n';
  T i = 0;
  for(auto a: inputs.pts) {
    std::cout << i <<" ("  << a.first << ","    << a.second
              << ")\t DEMAND:" << inputs.demand[i]  << "\n";
    i++;    //NOTE
  }
}


template<typename T1, typename T2>
float calRouteValue(const std::vector<T1> &aRoute, const std::vector< std::vector<T2>> &distance, T1 depot=1){ //return always a float
  float routeVal = 0.0f;
  T1 prevPoint = 0; //First point in a route is depot

  for(auto aPoint: aRoute) {
    routeVal += distance[prevPoint][aPoint];
    prevPoint=aPoint;
  }
  routeVal += distance[prevPoint][0]; //Last point in a route is depot

  return routeVal;
}


template<typename T1, typename T2>
void CheckPrintOutput(const CVRPInOut<T1> &inputs,
  const std::vector< std::vector<T2>> &distance
  ){

  float solVal = 0.0f;
  T1 i = 0;

  for(auto & aRoute: inputs.routesVector) { //NOTE

    float routeVal = 0.0f;
    T1 prevPoint = 0; //First point in a route is depot

    DEBUG std::cout<< "Route #"<< i++ <<":";  //increments only if prints

    for(auto aPoint: aRoute) {
      DEBUG std::cout<< " "<< aPoint ;
      routeVal += distance[prevPoint][aPoint];
      prevPoint=aPoint;
    }

    routeVal += distance[prevPoint][0]; //Last point in a route is depot

    DEBUG std::cout<< '\n';
    solVal += routeVal; //  can be done inside above for-block, but as of now, for clarity and it hardly matters as output is sparse

  }

  DEBUG std::cout<< "CostinOutput "<< inputs.cost << '\n'; //May need to round
  DEBUG std::cout<< "CostComputed "<< solVal      << '\n'; //May need to round
  if(fabs(inputs.cost-solVal) < 1.0f)
    std::cout<<  "MATCH"    << '\n';
  else {
    std::cout<<  "WARN: NOMATCH "  << fabs(inputs.cost-solVal)<< '\n';
    std::cout<< inputs.cost << " != "<< solVal << '\n';
  }
}

class Edge {
public:
	int to;
	float length;

	Edge(){}
	~Edge(){}
	Edge(int t, float l){
		to = t; length = l;
	}
	bool operator < (const Edge& e){
		return length < e.length;
	}
};

// This not templated as of now
// Modified W map to Dist[][] 2D
// CLRS version + tinkered
template <typename T1, typename T2>
std::vector<std::vector<Edge>>
PrimsAlgo( const std::vector<std::vector<Edge>> & graph, const std::vector< std::vector<T2>> &dist, T1 src){ //std::map< std::pair <int,int> , int> W
	auto N = graph.size();
  const T1 INIT = -1;
  //! std::cout<< "N "<< N << '\n';

	std::vector <T2> key(N, INT_MAX);
	std::vector <T2> toEdges(N, -1 );
	std::vector <bool> visited(N, false);

	std::set< std::pair<T2, T1> > active; // holds value and vertex


	//! key[0] = INT_MAX;
	//! visited[0] = true; // incorrect to set here!

	key[src] = 0.0;
	active.insert( { 0.0, src});

	while( active.size() > 0 ){
		auto where = active.begin()->second;

		//! DEBUG std::cout << "picked " << where <<"\tsize"<< active.size()<< std::endl;
		active.erase(active.begin());
		if(visited[where]) {
			continue;
		}
		visited[where] = true;
		for(Edge E : graph[where]){
			if(!visited[E.to] && E.length < key[E.to]){ //W[{where,E.to}]
				key[E.to] = E.length; //W[{where,E.to}]
				active.insert( { key[E.to], E.to});
				//! DEBUG std::cout << key[E.to] <<" ~ " <<  E.to << std::endl;
				toEdges[E.to]=where;
			}
		}
	}

	std::vector<std::vector<Edge>> nG(N);
	T1 u=0;
	for(auto v : toEdges){ // nice parallel code or made to parallel
		if(v != INIT ){
			//! int w = W[{u,v}];
			T2 w = dist[u][v];
			nG[u].push_back(Edge(v,w));
			nG[v].push_back(Edge(u,w));
			//! edges.push_back(std::make_pair(u,v));
		}
		u++;
	}
	return nG;
}

template <typename T1, typename T2>
std::vector<std::vector<Edge>> constructGraph(const CVRPInOut<T1> &inputs, const std::vector< std::vector<T2>> &dist){
  std::vector<std::vector<Edge>> nG(inputs.n);

  for(auto u = 0, endU = inputs.n; u < endU; ++u){
    for(auto v = 0, endV = inputs.n; v < endV; ++v){
      if(u < v){
        T2 w = dist[u][v];
        nG[u].push_back(Edge(v,w));
        nG[v].push_back(Edge(u,w));
      }
    }
  }

  return nG;
}

void printAdjList(const std::vector< std::vector<Edge> > &graph){
	int i = 0;
	for (auto vec : graph){

		std::cout << i << ": ";
		for(auto e : vec){
			std::cout<< e.to << " ";
		}
		i++;
		std::cout << std::endl;
	}
}

template <typename T, typename T2>
void ShortCircutTour(std::vector< std::vector<Edge> > &g, std::vector <bool> &visited, T u, std::vector<T2> &out){
  visited [u] = true;
  DEBUG std::cout<< u << ' ';
  //! cvrpio.addRouteVertex(u);
  out.push_back(u);
  for( auto e: g[u]){
    T v = e.to;
    if(!visited [v]){

      ShortCircutTour(g,visited,v, out);
    }
  }
}

template <typename T1, typename T2>
void convertToVrpRoutes(CVRPInOut<T1> &io, const std::vector <T2> &singleRoute) {

  T2 vCapacity  = io.vCapacity;
  T2 residueCap = vCapacity;
  std::vector<T1> aRoute;

  for (auto v : singleRoute){
    if (v == 0) continue;
    if(residueCap - io.demand[v] >= 0 ){
      aRoute.push_back(v);
      residueCap = residueCap - io.demand[v];
    }
    else{   //new route
      io.routes.push_back(aRoute);
      aRoute.clear();
      aRoute.push_back(v);
      residueCap = vCapacity - io.demand[v];
    }
  }
  io.routes.push_back(aRoute); //final route
}

template <typename T1, typename T2>
std::vector<std::vector<Edge>>
Dijkstra(const std::vector< std::vector<Edge> > &graph,
  const std::vector< std::vector<T2>> &dist,
  std::vector <T1> &singleRoute,
  T1 source = 0
	) {

  auto  N = graph.size();
  const T1 INIT = -1;

  std::vector <T1> parent(N, INIT);
	std::vector <T2> min_distance(N, (INT_MAX * 1.0)/3.0);

	min_distance[ source ] = 0.0;
	std::set< std::pair<T2,T1> > active_vertices; // (dist, vId)
	active_vertices.insert( {0.0,source} );

	while (!active_vertices.empty()) {
		T1 where = active_vertices.begin()->second;
    //! std::cout<< "(uv):" << where << " wt:" << active_vertices.begin()->first << '\n';
    singleRoute.push_back(where);
    active_vertices.erase( active_vertices.begin() );

		for (auto ed : graph[where]) {
			//! auto newdist = min_distance[where] + ed.length;
			auto newdistModified = min_distance[where] + ed.length + dist[0][where];
			if (newdistModified < min_distance[ed.to]) {
				active_vertices.erase( { min_distance[ed.to], ed.to } );
				min_distance[ed.to] = newdistModified;
				parent[ed.to] = where;
				active_vertices.insert( { newdistModified, ed.to } );
			}
		}
	}

  // Not necessary
  std::vector<std::vector<Edge>> nG(N);
	//! T1 u=0;
	//! for(auto v : parent){ // nice parallel code or made to parallel
		//! if(v != INIT ){
			//! T2 w = dist[u][v];

      //std::cout<< "(uv):" << u << " "<< v << " :"<< w << '\n';

			//! nG[u].push_back(Edge(v,w));
			//! nG[v].push_back(Edge(u,w));
		//! }
		//! u++;
	//! }
	return nG;
}



template<typename T>
void readVRPOutput(CVRPInOut<T> &inputs){
  std::string line, code, dummy;
  //! T depot = DEFAULT;
  T val = DEFAULT;
  //std::vector<T> allNodes;
  inputs.allNodes.push_back(0);
  while (std::getline(std::cin, line)) {
    std::istringstream iss(line);
    //! std::cout<< line << '\n';
    if(iss >> code && code != ""){
      //! std::cout<< code << '\n';
      transform(code.begin(), code.end(), code.begin(), ::toupper); // just in case if not upper
      if(code == "TIME"){
        // do nothing
      }else if (code == "ROUTE") {
        iss >> dummy;
        // std::cout<< "dum:" << dummy << '\n';
        std::vector<T> aRoute;
        while(iss >> val){
          //~ std::cout<< " "<< val ;
          aRoute.push_back(val);
          inputs.allNodes.push_back(val);
        }
        T sumDemand = 0;
        for(auto v: aRoute)
          sumDemand += inputs.demand[v];
        if(sumDemand > inputs.vCapacity)
          std::cout<< "WARN: Route demand sum > Capacity" << '\n';

        inputs.routesVector.push_back(aRoute);
        // std::cout<< " " << '\n';
      }else if (code == "COST"){
        iss >> inputs.cost;
      }else{
        std::cout<< "NEVER REACHED" << '\n';
      }
    }
  }
}
template<typename T>
void checkAll(CVRPInOut<T> &inputs){

  std::vector<T> allNodes ;
  std::vector<T> nodesRoute(inputs.n,DEFAULT);
  auto nGot   = inputs.allNodes.size();
  auto nGiven = inputs.n;
  if(nGot < (unsigned)nGiven)
    std::cout<< "WARN: Some nodes are missing " << nGot << " < " << nGiven << '\n';

  if(nGot > (unsigned)nGiven)
    std::cout<< "WARN: Some nodes are extraa " << nGot << " > " << nGiven << '\n';


  T routeId = 0;

  for(auto &aRoute: inputs.routesVector){
    DEBUG std::cout<< "R# "<< routeId + 1 << ": "; // +1 only for print purpose
    T aRouteCapacity = 0;
    for (auto &i: aRoute){
      aRouteCapacity += inputs.demand[i];
      DEBUG std::cout<< i<< " " ;
    }
    DEBUG std::cout<< "" << '\n';
    //check if each route is respecting vehicle capacity! - done in readOutput
    //~ if(abs(aRouteCapacity - inputs.vCapacity) < 1.0f){
      //~ std::cout << "NOT VALID: ROUTE CAP EXCEEDED for R#"<< routeId+1
                //~ << " AbsV:" << abs(aRouteCapacity - inputs.vCapacity) << '\n'
                //~ << " rCap:" << aRouteCapacity << '\n'
                //~ << " vCap:" << inputs.vCapacity << '\n';
    //~ }
    routeId++;
  }
}
int main(int argc, char **argv) {
  if(argc > 1){
    std::cout<< argv[1] << '\n';
    if (strcmp(argv[1],"-d")==0){
      DEBUGCODE = 1;
      std::cout<< "DEBUG MODE ON" << '\n';
    }
  }
  CVRPInOut  <float>  inputs;

  readVRP<float>(inputs); //input might be float! e.g. 33.000000



  std::vector< std::vector<float>> distance;
  calDist<float,float>(inputs,     distance);

  readVRPOutput<float>(inputs);   // Read the output file

  //~ std::cout<< "CAPTURED N:" << inputs.allNodes.size() << '\n';

  checkAll<float>(inputs);

  // Check if all vertices are present - Done
  // Check if there is no vertices present in two routes - Done
  // Check if each routes is respecting capacity - Done
  // Check if the cost output matches with the cost computed. - Done
  CheckPrintOutput<float>(inputs, distance);
  return 0;
}
