#ifndef _INTERACTION_NETWORKS_H_
#define _INTERACTION_NETWORKS_H_

#include <functional>
#include <map>
#include <limits>

#include "base/vector.h"
#include "base/assert.h"
#include "tools/Random.h"
#include "tools/random_utils.h"
#include "tools/vector_utils.h"
#include "tools/map_utils.h"
#include "tools/string_utils.h"
#include "tools/set_utils.h"
#include "tools/math.h"
#include "tools/distances.h"
#include "tools/attrs.h"
#include "tools/Graph.h"
#include "Evolve/Resource.h"

using org_t = emp::vector<int>;
using fit_map_t = std::map<org_t, double>;


DEFINE_ATTR(SigmaShare);
DEFINE_ATTR(Alpha);
DEFINE_ATTR(Cost);
DEFINE_ATTR(Cf);
DEFINE_ATTR(NicheWidth);
DEFINE_ATTR(MaxScore);
DEFINE_ATTR(ResourceInflow);
DEFINE_ATTR(ResourceOutflow);
DEFINE_ATTR(MaxBonus);

constexpr auto DEFAULT{MakeAttrs(SigmaShare(8.0),
                                 Alpha(1.0),
                                 Cost(1.0),
                                 Cf(.0025),                                 
                                 NicheWidth(3.0),
                                 MaxScore(10.0),
                                 ResourceInflow(2000),
                                 ResourceOutflow(.01),
                                 MaxBonus(5))};

using all_attrs = emp::tools::Attrs<typename SigmaShare::value_t<double>, typename Alpha::value_t<double>, 
                        typename Cost::value_t<double>, typename Cf::value_t<double>,
                        typename NicheWidth::value_t<double>, typename MaxScore::value_t<double>,
                        typename ResourceInflow::value_t<double>, typename ResourceOutflow::value_t<double>,
                        typename MaxBonus::value_t<double> >;


// from https://en.cppreference.com/w/cpp/types/numeric_limits/epsilon
template<class T, class U>
typename std::enable_if<!std::numeric_limits<T>::is_integer || !std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) <= std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

// integer version - just use the equals operator
template<class T, class U>
typename std::enable_if<std::numeric_limits<T>::is_integer && std::numeric_limits<U>::is_integer, bool>::type
    almost_equal(T x, U y, int ulp)
{
    return x == y;
}

/// Find the elements in @param pop with the highest value on
/// @param axis. PHEN_T should be a vector of some numeric type.
/// If epsilon lexicase is being used, specify a value of @param epsilon
/// to include everything within that range of the highest value
template <typename PHEN_T>
emp::vector<PHEN_T> FindHighest(emp::vector<PHEN_T> & pop, int axis, double epsilon = 0) {
    double best = -1;
    emp::vector<PHEN_T> winners;

    for (PHEN_T & org : pop) {
        if (org[axis] > best) {
            best = org[axis];
            winners.resize(0);
            winners.push_back(org);
        } else if (almost_equal(org[axis], best, 5)) {
            winners.push_back(org);
        }
    }

    if (epsilon) {
        winners.resize(0);
        for (PHEN_T & org : pop) {
            if (org[axis] + epsilon >= best) {
                winners.push_back(org);
            }
        }
    }
    return winners;
}

/// Multiply all of the ints in @param nums togetehr.
long int VectorProduct(emp::vector<int> nums) {
    long int result = 1;

    for (int num : nums) {
        result *= num;
    }

    return result;
}

template <typename PHEN_T>
emp::vector<PHEN_T> FilterImpossible(emp::vector<PHEN_T> & pop, emp::vector<int> & axes, double epsilon = 0) { 
    std::set<PHEN_T> possible;
    
    for (int ax : axes) {
        emp::vector<PHEN_T> best = FindHighest(pop, ax, epsilon);
        possible = emp::set_union(possible, best);
    }

    emp::vector<PHEN_T> result(possible.begin(), possible.end());
    return result;
}

template <typename PHEN_T>
emp::vector<int> PruneAxes(emp::vector<int> & axes, emp::vector<PHEN_T> & pop, double epsilon = 0) {
    emp::vector<int> result;

    for (int ax : axes) {
        // std::cout << "Evaluating " << ax << std::endl;
        bool include = false;
        double best = pop[0][ax];
        double lowest = pop[0][ax];
        for (PHEN_T & org : pop) {
            if (org[ax] > best) {
                // std::cout << org[ax] << " > "<< best << std::endl;
                if (org[ax] > lowest + epsilon) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                best = org[ax];
            } else if (org[ax] < lowest) {
                // std::cout << org[ax] << " < "<< lowest << std::endl;
                if (org[ax] + epsilon < best) {
                    // std::cout <<" Including  "<< ax << std::endl;
                    include = true;
                    break;
                }
                lowest = org[ax];
            }
        }

        if (include) {
            result.push_back(ax);
        }
    }

    return result;

}

template <typename PHEN_T>
void TraverseDecisionTree(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & pop, emp::vector<int> axes, emp::vector<int> perm_levels, double epsilon = 0) {
    // std::cout << "begining round of recursion " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    if (axes.size() == 1) {
        emp::vector<PHEN_T> winners = FindHighest(pop, axes[0], epsilon);
        // std::cout << "Winners: " << emp::to_string(winners) << std::endl;
        for (PHEN_T & org : winners) {
            fit_map[org]+=1.0/((double)winners.size()*VectorProduct(perm_levels));
        }
        return;
    }

    axes = PruneAxes(axes, pop, epsilon);
    perm_levels.push_back(axes.size());

    pop = FilterImpossible(pop, axes, epsilon);

    // std::cout << "post processing: " << axes.size() << emp::to_string(pop) << emp::to_string(axes) << std::endl;

    for (int ax : axes) {
        // std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        emp::vector<PHEN_T> winners = FindHighest(pop, ax, epsilon);
        if (winners.size() == 1) { // Not a tie
            // std::cout << "1 winner: " << emp::to_string(winners[0]) << " Controls " << (double)emp::Factorial(axes.size() - 1)<< std::endl;
            fit_map[winners[0]] += 1.0/VectorProduct(perm_levels);//(double)emp::Factorial(axes.size() - 1);
        } else { // tie
            emp::vector<int> next_axes;
            for (int new_ax : axes) {
                if (new_ax != ax) {
                    next_axes.push_back(new_ax);
                }
            }
            TraverseDecisionTree(fit_map, winners, next_axes, perm_levels, epsilon);
        }
    }
}


template <typename PHEN_T>
std::map<PHEN_T, double> LexicaseFitness(emp::vector<PHEN_T> & pop, double epsilon = 0) {

    emp_assert(pop.size() > 0);
    std::map<PHEN_T, double> fit_map;
    size_t n_funs = pop[0].size();

    for (PHEN_T & org : pop) {
        fit_map[org] = 0.0;
    }

    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTree(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs), {}, epsilon);

    for (PHEN_T & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
        // std::cout << "Pre div: " << fit_map[org] << std::endl;
        // fit_map[org] /= emp::Factorial(n_funs); // convert to proportion of "islands"
        // std::cout << "Post div: " << fit_map[org] << " Divided by: " << (double)emp::Factorial(n_funs) << std::endl;
    }

    return fit_map;
}

template <typename PHEN_T>
void TournamentHelper(std::map<PHEN_T, double> & fit_map, emp::vector<PHEN_T> & pop, int t_size = 2){

    std::map<PHEN_T, double> base_fit_map = fit_map;

    for (org_t & org : pop) {
        double less = 0.0;
        double equal = 0.0; 

        for (org_t & org2 : pop) {
            if (almost_equal(base_fit_map[org2], base_fit_map[org], 10)) {
                equal++;
            } else if (base_fit_map[org2] < base_fit_map[org]) {
                less++;
            }
        }

        long double p_less = less/(double)pop.size();
        long double p_equal = equal/(double)pop.size();
        long double p_self = 1/(double)pop.size();

        fit_map[org] = (pow(p_equal + p_less, t_size) - pow(p_less, t_size)) * p_self/p_equal;
    }
}

template <typename PHEN_T>
std::map<PHEN_T, double> TournamentFitness(emp::vector<PHEN_T> & pop, int t_size = 2) {
    emp_assert(pop.size() > 0);
    std::map<PHEN_T, double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map[org] = emp::Sum(org);
    }
    TournamentHelper(fit_map, pop, t_size);
    return fit_map;
}


std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> sharing_fitness = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "SHARING" << std::endl;
    fit_map_t fit_map;

    // std::cout << SigmaShare::Get(attrs) << std::endl;

    for (org_t & org : pop) {
        fit_map[org] = 1.0;
    }

    for (org_t & org1 : pop) {
        double niche_count = 0;

        for (org_t & org2 : pop) {
            // Sharing function is euclidean distance
            // we could make this configurable
            double dist = emp::EuclideanDistance(org1, org2);
            if (dist < SigmaShare::Get(attrs)) {
                niche_count += 1 - emp::Pow((dist/SigmaShare::Get(attrs)), Alpha::Get(attrs));
            } 
        }

        // Don't worry, niche_count will always be at least one because each individual
        // increases its own niche count by 1
        fit_map[org1] = emp::Sum(org1) / niche_count;
    }
    TournamentHelper(fit_map, pop, 2);

    return fit_map;    
};


std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> eco_ea_fitness = [](emp::vector<org_t> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "ECO-EA" << std::endl;
    emp_assert(pop.size() > 0);

    fit_map_t base_fit_map;
    fit_map_t fit_map;

    size_t n_funs = pop[0].size();

    for (org_t & org : pop) {
        fit_map[org] = 1.0;
    }

    for (int axis : emp::NRange(0, (int)n_funs)) {
        double res = ResourceInflow::Get(attrs);
        double count = 0;

        for (org_t & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                count++;
            }
        }

        if (count > 0) {
            res /= count; // Ignores resource accumulation, but that's okay for interactio  }
        }
        for (org_t & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                fit_map[org] *= emp::Pow2(std::min(Cf::Get(attrs)*res*emp::Pow(org[axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs)));
            }   
        }
    }

    TournamentHelper(fit_map, pop, 2);

    return fit_map;
};

emp::WeightedGraph calc_competition(emp::vector<org_t> pop, 
                        std::function<fit_map_t(emp::vector<org_t>&, all_attrs)> fit_fun,
                        all_attrs attrs=DEFAULT) {

    emp::WeightedGraph effects(pop.size());


    fit_map_t fitnesses = fit_fun(pop, attrs);

    for (size_t i = 0; i < pop.size(); i++) {
        effects.SetLabel(i, emp::to_string(pop[i]));
        // std::cout << effects.GetLabel(i) << std::endl;

        emp::vector<org_t> curr = pop;
        for (int & ax : curr[i]) {
            ax = 0; // Replace org with null org so pop size stays same
        }

        fit_map_t new_fits = fit_fun(curr, attrs);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}
            double effect = fitnesses[curr[j]] - new_fits[curr[j]];
            if (effect != 0) {
                effects.AddEdge(i, j, effect);
            }
            // std::cout << effect << std::endl;
        }
    }

    return effects;
}


#endif
