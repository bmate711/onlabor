#include "ortools/base/logging.h"
#include "ortools/constraint_solver/routing.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>

namespace operations_research
{
class ReadData
{
  private:
    int nWP;
    int nConfig;
    int nJoint;
    std::vector<std::vector<std::vector<double> *> *> *graspConfigs;
    std::vector<std::vector<double> *> *putawayConfig;
    std::vector<double> *endConfig;
    std::vector<double> *startConfig;
    
    std::vector<double> *split(std::string line)
    {
        std::vector<double> *node = new std::vector<double>();
        int pos;
        while ((pos = line.find("\t")) != std::string::npos)
        {
            double angle = atof(line.substr(0, pos).c_str());
            node->push_back(angle);
            line.erase(0, pos + 1);
        }
        double angle = atof(line.substr(0, pos).c_str());
        node->push_back(angle);
        return node;
    }

  public:
    ReadData(int i)
    {
        std::string fileName = "D:\\Prog\\C++\\or-tools\\onlab\\Instances\\test";
        fileName += std::to_string(i);
        fileName += ".txt";
        std::string line;
        std::ifstream data(fileName);

        if (data.is_open())
        {
            getline(data, line);
            nWP = std::stoi(line);
            
            getline(data, line);
            nConfig = std::stoi(line);

            getline(data, line);
            nJoint = std::stoi(line);

            getline(data, line);
            startConfig = split(line);

            getline(data, line);
            endConfig = split(line);

            putawayConfig = new std::vector<std::vector<double> *>();
            for (int i = 0; i < nConfig; ++i)
            {
                getline(data, line);
                putawayConfig->push_back(split(line));
            }
            // for(int i = 0; i < nConfig-4; ++i){
            //     std::vector<double> *node = new std::vector<double>();
            //     for(int j = 0; j < nJoint; ++j) {
            //         node->push_back(9999.9);
            //     }
            //     putawayConfig->push_back(node);
            // }

            graspConfigs = new std::vector<std::vector<std::vector<double> *> *>();
            std::vector<std::vector<double> *> *group = new std::vector<std::vector<double> *>();
            int i = 1;
            while (getline(data, line))
            {
                group->push_back(split(line));
                if (i % nConfig == 0)
                {
                    graspConfigs->push_back(group);
                    group = new std::vector<std::vector<double> *>();
                }
                i++;
            }
        }
    }

    int GetNWP() const
    {
        return nWP;
    }

    int GetNConfig() const
    {
        return nConfig;
    }

    int GetNJoint() const
    {
        return nJoint;
    }

    std::vector<double> *GetEndConfig() const
    {
        return endConfig;
    }

    std::vector<double> *GetStartConfig() const
    {
        return startConfig;
    }

    std::vector<double> &GetConfig(int index) const
    {
        if (index == 0)
            return *startConfig;
        if (index < nWP * nConfig + 1)
            return *graspConfigs->at((index - 1) / nConfig)->at((index - 1) % nConfig);
        // if (index > nWP * nConfig)
        return *putawayConfig->at((index - 1) % nConfig);
    }

    std::vector<std::vector<std::vector<double> *> *> *GetGraspConfig() const
    {
        return graspConfigs;
    }

    std::vector<std::vector<double> *> *GetPutawayConfig() const
    {
        return putawayConfig;
    }

    std::size_t GetVehicleNumber() const { return 1; }

    RoutingModel::NodeIndex GetDepot() const { return RoutingModel::kFirstNode; }

    ~ReadData()
    {
        for (auto group : *graspConfigs)
        {
            for (auto node : *group)
                delete node;
            delete group;
        }
        delete graspConfigs;

        for (auto node : *putawayConfig)
            delete node;
        delete putawayConfig;

        delete startConfig;

        delete endConfig;
    }
};

class RotationTimeMatrix : public RoutingModel::NodeEvaluator2
{
  private:
    std::vector<double> acceleratios;
    std::vector<double> speed;
    std::vector<std::vector<int64>> distances;
    int64 infiniteWeight = 999999999999999;
    ReadData *readData;
    std::vector<std::vector<RoutingModel::NodeIndex> *> *disjunctions;

    void CreateDisjunction()
    {
        disjunctions = new std::vector<std::vector<RoutingModel::NodeIndex> *>();
        for (int i = 0; i < readData->GetNWP() * 2; ++i)
        {
            auto nodes = new std::vector<RoutingModel::NodeIndex>();
            for (int j = 0; j < readData->GetNConfig(); ++j)
            {
                nodes->push_back((RoutingModel::NodeIndex)(i * readData->GetNConfig() + j + 1));
            }
            disjunctions->push_back(nodes);
        }
    }

    int64 CalculateWeight(std::vector<double> *from, std::vector<double> *to)
    {
        double max = 0.0;
        for (int i = 0; i < readData->GetNJoint(); ++i)
        {

            double angleDiff = abs(from->at(i) - to->at(i));
            double time = 0;
            if (IsTriangleProfile(angleDiff, speed.at(i), acceleratios.at(i)))
            {
                time = sqrt(4 * angleDiff / acceleratios.at(i));
            }
            else
            {
                time = TrapezoidProfile(angleDiff, speed.at(i), acceleratios.at(i));
            }
            if (max < time)
                max = time;
        }
        return Double2Int64(max);
    }

    int64 Double2Int64(double number)
    {
        return (int64)((number)*10000000 + 0.5); //10^7
    }

    bool IsTriangleProfile(double angleDiff, double speed, double acceleration)
    {
        return angleDiff < speed * speed / acceleration;
    }

    double TrapezoidProfile(double angleDiff, double speed, double acceleration)
    {
        return 2.0 * speed / acceleration + (angleDiff - (speed * speed / acceleration)) / speed;
    }

    void CreateMatrix()
    {
        StartWeigth();
        GrepsWeigth();
        PutawayWeight();
        // EndWeight();
    }

    // void EndWeight(){
    //     auto row = std::vector<int64>();
    //     for(int i = 0; i < 2 * (readData->GetNWP() * readData->GetNConfig() + 1); ++i){
    //         row.push_back(infiniteWeight);
    //     }
    //     distances.push_back(row);
    // }

    void PutawayWeight()
    {
        for (int i = 0; i < readData->GetNWP(); ++i)
        {
            for (auto from : *readData->GetPutawayConfig())
            {
                auto row = std::vector<int64>();
                row.push_back(CalculateWeight(from, readData->GetStartConfig())); //row.push_back(infinteWeight)
                for (auto group : *readData->GetGraspConfig())
                {
                    for (auto to : *group)
                    {
                        row.push_back(CalculateWeight(from, to));
                    }
                }
                for (int j = 0; j < readData->GetNWP() * readData->GetNConfig(); ++j)
                {
                    row.push_back(infiniteWeight);
                }
                // row.push_back(CalculateWeight(from, readData->GetEndConfig()));
                distances.push_back(row);
            }
        }
    }

    void StartWeigth()
    {
        auto row = std::vector<int64>();
        row.push_back(infiniteWeight);
        for (auto group : *readData->GetGraspConfig())
        {
            for (auto node : *group)
            {
                row.push_back(CalculateWeight(readData->GetStartConfig(), node));
            }
        }
        for (int i = 0; i < readData->GetNWP() * readData->GetNConfig(); ++i)
        { //+1
            row.push_back(infiniteWeight);
        }
        this->distances.push_back(row);
    }

    void GrepsWeigth()
    {

        for (auto group : *readData->GetGraspConfig())
        {
            for (auto from : *group)
            {
                auto row = std::vector<int64>();
                for (int i = 0; i < readData->GetNWP() * readData->GetNConfig() + 1; ++i)
                {
                    row.push_back(infiniteWeight);
                }
                for (int i = 0; i < readData->GetNWP(); ++i)
                {
                    for (auto to : *readData->GetPutawayConfig())
                    {
                        row.push_back(CalculateWeight(from, to));
                    }
                }
                // row.push_back(infiniteWeight);
                distances.push_back(row);
            }
        }
    }

  public:
    RotationTimeMatrix(ReadData *readData)
    {
        this->readData = readData;
        for (int i = 0; i < readData->GetNJoint(); ++i) //TODO: Read from file
        {
            acceleratios.push_back(10.472);
            if (i < 3)
                speed.push_back(2.269);
            else
                speed.push_back(3.316);
        }
        CreateMatrix();
        CreateDisjunction();
    }

    std::vector<std::vector<RoutingModel::NodeIndex> *> *GetDisjuntions()
    {
        return disjunctions;
    }

    int64 Run(RoutingModel::NodeIndex FromNode,
              RoutingModel::NodeIndex ToNode) override
    {
        return distances[FromNode.value()][ToNode.value()];
    }

    ~RotationTimeMatrix()
    {
        for(auto i : *disjunctions){
            delete i;
        }
        delete disjunctions;
    }
};

void PrintSolution(const ReadData &data, const RoutingModel &routing,
                   const Assignment &solution)
{
    int64 index = routing.Start(0);  
    int64 distance = 0LL;
    std::stringstream route;
    int nodeIndex;
    auto allIndex = std::vector<int>();
    int q = 1;
    while (routing.IsEnd(index) == false)
    {
        nodeIndex = routing.IndexToNode(index).value();
        route << nodeIndex << " " << (nodeIndex - 1 )/data.GetNConfig() + 1 << " " << nodeIndex % data.GetNConfig();
        for (auto i : data.GetConfig(nodeIndex))
        {
            route << " " << i;
        }
        route << std::endl;
        int64 previous_index = index;
        index = solution.Value(routing.NextVar(index));
        //if(q != 1 && q != data.GetNWP()*2+1)
        // std::cout << (double)const_cast<RoutingModel &>(routing).GetArcCostForVehicle(
        //          previous_index, index, 0LL)/10000000.0 << std::endl; 
        distance += const_cast<RoutingModel &>(routing).GetArcCostForVehicle(
                previous_index, index, 0LL);
        q++;
    }

    nodeIndex = routing.IndexToNode(index).value();
    route << nodeIndex << " " << (nodeIndex - 1 )/data.GetNConfig() + 1 << " " << nodeIndex % data.GetNConfig();
    for (auto i : data.GetConfig(nodeIndex))
    {
        route << " " << i;
    }
    route<<std::endl;
    // std::cout << route.str();
    // std::cout << "Distance of the route: " << ((double)distance / 10000000) << "s" << std::endl;
    // std::cout << "Advanced usage:";
    // std::cout << "Problem solved in " << routing.solver()->wall_time() << "ms";
    // std::cout << std::endl << std::endl;
    //std::cout << routing.solver()->wall_time() << std::endl;
    std::cout<< ((double)distance / 10000000)<<std::endl;
}

void Solve(int i)
{
    try{
    ReadData *readData = new ReadData(i);
    RoutingModel routing(readData->GetNWP() * readData->GetNConfig() * 2 + 1, readData->GetVehicleNumber(), readData->GetDepot());
    RotationTimeMatrix *matrix = new RotationTimeMatrix(readData);
    for (auto row : *matrix->GetDisjuntions())
    {
        routing.AddDisjunction(*row);
    }
    routing.SetArcCostEvaluatorOfAllVehicles(NewPermanentCallback(matrix, &RotationTimeMatrix::Run));
    RoutingSearchParameters searchParameters = RoutingModel::DefaultSearchParameters();
    searchParameters.set_first_solution_strategy(FirstSolutionStrategy::PATH_CHEAPEST_ARC);
    searchParameters.set_local_search_metaheuristic(LocalSearchMetaheuristic::TABU_SEARCH);
    searchParameters.set_time_limit_ms(20);
    const Assignment *solution = routing.SolveWithParameters(searchParameters);
    // const Assignment* solution = routing.Solve();
    PrintSolution(*readData, routing, *solution);
    int64 index = routing.Start(0);
    std::stringstream route;
    delete matrix;
    delete readData;
    } 
    catch (std::exception& e){
        std::cout << e.what() << '\n';
    }
}
} // namespace operations_research

int main(int argc, char **argv)
{
    for(int i = 1; i <= 40; ++i ){
        operations_research::Solve(i);
        // std::cout<<std::endl<<std::endl;    
    }
    // operations_research::Solve(17);
    // operations_research::Solve(19);
}