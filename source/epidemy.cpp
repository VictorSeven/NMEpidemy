#include<iostream>
#include<cstdlib>
#include<cmath>
#include<fstream>
#include<random>
#include<vector>
#include<tuple>
#include<string>
#include "CNet.cpp"

#define MAX 1000

using namespace std;

void init(CNetwork net, bool infected[MAX], vector<int>& nodes, vector<double>& nodes_time);

tuple<int, double, double> gillespie(int Nrates, double rate_sum, vector<double> rates, mt19937& gen, uniform_real_distribution<double>& random);

double get_rates(double tau, vector<double>& nodes_time, vector<double>& links_time, vector<double>& result, double beta0, double lambda0);

void work_list(int mu, bool infected[MAX], CNetwork net, vector<int>& nodes, vector<double>& nodes_time, vector<int>& links, vector<double>& links_time);

double beta(double t, double beta0);
double lambda(double t, double lambda0);

int main(void)
{
    int i,j,k;

    double rate_sum;
    double beta0 = 1.0;
    double lambda0 = 0.4;

    int its = 3e4;

    mt19937 gen(958433198);
    uniform_real_distribution<double> ran_u(0.0,1.0);

    ofstream output;

    CNetwork net(MAX, false);



    int netsize;

    /*vector<int> nodes(netsize);
    vector<double> nodes_time(netsize);

    vector<int> links = vector<int>();
    vector<double> links_time = vector<double>();
    vector<double> rates = vector<double>();*/

    vector<int> nodes = vector<int>();
    vector<double> nodes_time = vector<double>();
    vector<int> links = vector<int>();
    vector<double> links_time = vector<double>();
    vector<double> rates = vector<double>();
    bool infected[MAX];

    double mean_pop;
    int n_measures;

    double avg_infected, sqr_infected;
    int averages = 50;

    double t, tau;
    double a;
    int mu;

    output.open("constant.txt");

    for (lambda0 = 0.00; lambda0 < 1; lambda0 += 0.02)
    {
        avg_infected = sqr_infected = 0.0;
        for (j=0; j < averages; j++)
        {
            net = CNetwork(MAX, false);
            net.create_erdos_renyi(MAX, 8.0, 958431190 + 10000 * j);
            netsize = net.get_node_count();

            //Init all the vector stuff
            nodes = vector<int>(netsize);
            nodes_time = vector<double> (netsize);
            links = vector<int>();
            links_time = vector<double>();
            rates = vector<double>();

            init(net, infected, nodes, nodes_time); //Start algorithm

            tau = 0.0; //Init tau to 0

            //Init time and also things for average
            t = 0.0;
            mean_pop = 0.0;
            n_measures = 0;
            for (i=0; i < its; i++)
            {
                rate_sum = get_rates(tau, nodes_time, links_time, rates, beta0, lambda0); //Get total and partial rates
                tie(mu, tau, a) = gillespie(netsize, rate_sum, rates, gen, ran_u); //Execute gillespie
                t += tau; //Update time
                //Make changes in list
                work_list(mu, infected, net, nodes, nodes_time, links, links_time);
                //output << t << " " << nodes.size() << endl;

                //For enough time passed, make averages
                if (t >= 10)
                {
                    mean_pop += nodes.size(); //Number of infected
                    n_measures += 1;
                }

                //If there is no people, go out of loop
                if (nodes.size() == 0)
                {
                    mean_pop = 0.0;
                    n_measures = 1;
                    break;
                }
            } //End i for

            mean_pop /= (1.0*n_measures);
            avg_infected += mean_pop;
            sqr_infected += mean_pop * mean_pop;

        } //End j for

    avg_infected /= 1.0 * averages;
    sqr_infected /= 1.0 * averages;

    //Finish and write averages

    output << lambda0/beta0 << " " << avg_infected << " " << sqrt((sqr_infected - avg_infected* avg_infected)/(1.0*averages)) << endl;

    }

    output.close();


    return 0;
}

//Init the thing
void init(CNetwork net, bool infected[MAX], vector<int>& nodes, vector<double>& nodes_time)
{
    int i;
    vector<int> who_to_erase = vector<int>();
    for (i=0; i < net.get_node_count(); i++)
    {
        //Do I have a disconnected node?
        if (net.degree(i) == 0)
        {
            who_to_erase.push_back(i); //Take note to erase

            infected[i] = false; //Susceptible forever
        }
        else //In other case, all nodes infected
        {
            nodes[i] = i;
            nodes_time[i] = 0.0;
            infected[i] = true;
        }

    }


    //Make the deletion of the nodes, taking in account that this changes indices
    for(i=0; i < who_to_erase.size(); i++)
    {
        nodes.erase(nodes.begin() + (who_to_erase[i]-i));
        nodes_time.erase(nodes_time.begin() + (who_to_erase[i]-i));
    }
}


double beta(double t, double beta0)
{
    double alpha = 2.0;
    //
    //return beta0 * exp(beta0 * t);
    return beta0;
}

double lambda(double t, double lambda0)
{
    double alpha = 0.5;
    //return (2*t)/(alpha + t*t)*sqrt(alpha)*lambda0;
    //return t == 0.0 ? 0.0 : alpha * pow(lambda0, alpha) * pow(t, alpha - 1);
    return lambda0;
}

double get_rates(double tau, vector<double>& nodes_time, vector<double>& links_time, vector<double>& result, double beta0, double lambda0)
{
    int i;
    double sum;

    //Make result able to store all the rates
    result = vector<double>(nodes_time.size() + links_time.size());

    //Compute all the rates from infected nodes,
    for (i=0; i < nodes_time.size(); i++)
    {
        sum += beta(nodes_time[i], beta0);
        nodes_time[i] += tau; //Also advance time
        result[i] = sum;
    }

    //and for infection links
    for (i=0; i < links_time.size(); i++)
    {
        sum += lambda(links_time[i], lambda0);
        links_time[i] += tau;
        result[nodes_time.size()+i] = sum;
    }

    return sum; //Return total rate
}

void work_list(int mu, bool infected[MAX], CNetwork net, vector<int>& nodes, vector<double>& nodes_time, vector<int>& links, vector<double>& links_time)
{
    int i,j,k;
    int index;
    int aux;


    vector<int> who_to_to_erase = vector<int>(); //Auxiliary list to erase well

    if (mu < nodes.size()) //Selected an infected node
    {
        index = nodes[mu];

        //Eliminate from list of infected
        nodes.erase(nodes.begin() + mu);
        nodes_time.erase(nodes_time.begin() + mu);

        infected[index] = false; //Make infected

        //Eliminate links between i and other
        for (i=0; i < links_time.size(); i++)
        {
            if (links[2*i] == index)
            {
                who_to_to_erase.push_back(2*i);
            }
        }

        //Erase all the things in the list
        for (i=0; i < who_to_to_erase.size(); i++)
        {
            links.erase(links.begin() + (who_to_to_erase[i] - 2*i));
            links.erase(links.begin() + (who_to_to_erase[i] - 2*i));
            links_time.erase(links_time.begin() + (who_to_to_erase[i]/2 - i));
        }

        //Add to the list connections between i and infected nodes
        for (i=0; i < net.get_num_neighs(index); i++)
        {
            aux = net.get_neigh_at(index, i);
            if (infected[aux])
            {
                links.push_back(aux);
                links.push_back(index);
                links_time.push_back( nodes_time[aux] ); //Rule 2: assign the time node has been infected
            }
        }
    }
    else
    {

        index = links[2 * (mu-nodes.size()) + 1];//Take the susceptible. They are always ordered as I->S

        infected[index] = true; //Infect

        //Add node to the list of infected
        nodes.push_back(index);
        nodes_time.push_back(0.0);

        //Eliminate links between i and other infected
        for (i=0; i < links_time.size(); i++)
        {
            if (links[2*i+1] == index)//Nodes is to. Then eliminate i-1 and i
            {
                who_to_to_erase.push_back(2*i);
            }

        }

        //Finally erase everything
        for (i=0; i < who_to_to_erase.size(); i++)
        {
            links.erase(links.begin() + (who_to_to_erase[i] - 2*i));
            links.erase(links.begin() + (who_to_to_erase[i] - 2*i));
            links_time.erase(links_time.begin() + (who_to_to_erase[i]/2 - i));
        }

        //Add to the list connections between i and infected nodes
        for (i=0; i < net.get_num_neighs(index); i++)
        {
            aux = net.get_neigh_at(index, i);
            if (not infected[aux])
            {
                links.push_back(index);
                links.push_back(aux);
                links_time.push_back(0.0);
            }
        }
    }
    return;
}

tuple<int, double, double> gillespie(int Nrates, double rate_sum, vector<double> rates, mt19937& gen, uniform_real_distribution<double>& random)
{
    int i; //Counter
    int mu; //Reaction to do

    double left, right; //To compute the inequalities

    double tau,a; //Parameters to return

    double r1, r2; //Two random numbers

    a = rate_sum; //Get the total rate

    //Generate the two random numbers
    r1 = random(gen);
    r2 = random(gen);

    tau =  -log(r1) / a; //Time

    //To compute the reaction...
    mu = 0;
    left = 0.0;
    right = rates[0] / a;

    //Get the mu using a binary search to locate which is the first partial rate greater than it
    //vector<double> partial_rates = vector<double>(rates, rates + Nrates);
    mu = distance(rates.begin(), upper_bound(rates.begin(), rates.end(), a*r2));


    //Get between random rate_parcial[mu] and rate_parcial[mu+1]
    /*while( not (left < r2 and r2 < right) )
    {
        mu += 1; //Increase counter
        left = right; //Update left
        right = rates[mu] / a; //Next to the right
    }*/

    return make_tuple(mu,tau,a); //Make the tuple and return it
}
