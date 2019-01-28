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

#include <math.h>  

DEFINE_ATTR(SigmaShare);
DEFINE_ATTR(Alpha);
DEFINE_ATTR(Cost);
DEFINE_ATTR(Cf);
DEFINE_ATTR(NicheWidth);
DEFINE_ATTR(MaxScore);
DEFINE_ATTR(ResourceInflow);
DEFINE_ATTR(ResourceOutflow);
DEFINE_ATTR(MaxBonus);
DEFINE_ATTR(TournamentSize);
DEFINE_ATTR(Epsilon);

constexpr auto DEFAULT{MakeAttrs(SigmaShare(8.0),
                                 Alpha(1.0),
                                 Cost(1.0),
                                 Cf(.0025),                                 
                                 NicheWidth(3.0),
                                 MaxScore(10.0),
                                 ResourceInflow(2000.0),
                                 ResourceOutflow(.01),
                                 MaxBonus(5.0),
                                 TournamentSize(2),
                                 Epsilon(0.0))};

using all_attrs = emp::tools::Attrs<typename SigmaShare::value_t<double>, typename Alpha::value_t<double>, 
                        typename Cost::value_t<double>, typename Cf::value_t<double>,
                        typename NicheWidth::value_t<double>, typename MaxScore::value_t<double>,
                        typename ResourceInflow::value_t<double>, typename ResourceOutflow::value_t<double>,
                        typename MaxBonus::value_t<double>, typename TournamentSize::value_t<int>, typename Epsilon::value_t<double> >;


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

/// Multiply all of the ints in @param nums together.
long int VectorProduct(const emp::vector<int> & nums) {
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
        for (PHEN_T & p : best) {
            possible.insert(p);
        }
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
        // if (perm_levels.size() == 0) {
        //     std::cout << "Axis: " << ax << " out of " << emp::to_string(axes) << std::endl;
        // }
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
emp::vector<double> LexicaseFitness(emp::vector<PHEN_T> & pop, all_attrs attrs = DEFAULT) {

    emp_assert(pop.size() > 0);
    std::map<PHEN_T, double> fit_map;
    size_t n_funs = pop[0].size();

    for (PHEN_T & org : pop) {
        fit_map[org] = 0.0;
    }

    emp::vector<PHEN_T> de_dup_pop = emp::RemoveDuplicates(pop);
    TraverseDecisionTree(fit_map, de_dup_pop, emp::NRange(0, (int)n_funs), {}, Epsilon::Get(attrs));

    for (PHEN_T & org : de_dup_pop) {
        fit_map[org] /= emp::Count(pop, org);
        // std::cout << "Pre div: " << fit_map[org] << std::endl;
        // fit_map[org] /= emp::Factorial(n_funs); // convert to proportion of "islands"
        // std::cout << "Post div: " << fit_map[org] << " Divided by: " << (double)emp::Factorial(n_funs) << std::endl;
    }

    emp::vector<double> result;
    for (PHEN_T & org : pop) {
        result.push_back(fit_map[org]);
    }

    return result;
}

void TournamentHelper(emp::vector<double> & fit_map, int t_size = 2){

    emp::vector<double> base_fit_map = fit_map;
    double pop_size = base_fit_map.size();

    for (int i = 0; i < fit_map.size(); i++) {
        double less = 0.0;
        double equal = 0.0; 

        for (auto & org2 : base_fit_map) {
            if (almost_equal(org2, base_fit_map[i], 10)) {
                equal++;
            } else if (org2 < base_fit_map[i]) {
                less++;
            }
        }

        long double p_less = less/pop_size;
        long double p_equal = equal/pop_size;
        long double p_self = 1/pop_size;

        fit_map[i] = (pow(p_equal + p_less, t_size) - pow(p_less, t_size)) * p_self/p_equal;
    }
}

template <typename PHEN_T>
emp::vector<double> TournamentFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    emp_assert(pop.size() > 0);
    // emp::vector<double> fit_map;
    emp::vector<double> fit_map;
    for (PHEN_T & org : pop) {
        fit_map.push_back(emp::Sum(org));
    }
    TournamentHelper(fit_map, TournamentSize::Get(attrs));
    return fit_map;
}

template <typename PHEN_T>
emp::vector<double> SharingFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "SHARING" << std::endl;
    emp::vector<double> fit_map;

    // std::cout << SigmaShare::Get(attrs) << std::endl;

    for (int i = 0; i < pop.size(); i++) {
        fit_map.push_back(1.0);
    }

    for (int i = 0; i < pop.size(); i++) {
        double niche_count = 0;

        for (PHEN_T & org2 : pop) {
            // Sharing function is euclidean distance
            // we could make this configurable
            double dist = emp::EuclideanDistance(pop[i], org2);
            if (dist < SigmaShare::Get(attrs)) {
                niche_count += 1 - pow((dist/SigmaShare::Get(attrs)), Alpha::Get(attrs));
            } 
        }

        // Don't worry, niche_count will always be at least one because each individual
        // increases its own niche count by 1
        fit_map[i] = emp::Sum(pop[i]) / niche_count;
    }
    TournamentHelper(fit_map, TournamentSize::Get(attrs));

    return fit_map;    
};

template <typename PHEN_T>
emp::vector<double> EcoEAFitness(emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT) {
    // std::cout << "ECO-EA" << std::endl;
    emp_assert(pop.size() > 0);

    emp::vector<double> fit_map;

    size_t n_funs = pop[0].size();

    for (int i = 0; i < pop.size(); i++) {
        fit_map.push_back(1.0);
    }

    for (int axis : emp::NRange(0, (int)n_funs)) {
        double res = ResourceInflow::Get(attrs);
        double count = 0;

        for (PHEN_T & org : pop) {
            if (org[axis] >= NicheWidth::Get(attrs)) {
                count++;
            }
        }

        if (count > 0) {
            res /= count; // Ignores resource accumulation, but that's okay for interactio  }
        }
        for (int i = 0; i < pop.size(); i++) {
            if (pop[i][axis] >= NicheWidth::Get(attrs)) {
                fit_map[i] *= pow(2,std::min(Cf::Get(attrs)*res*pow(pop[i][axis]/MaxScore::Get(attrs),2.0) - Cost::Get(attrs), MaxBonus::Get(attrs)));
            }   
        }
    }

    TournamentHelper(fit_map, TournamentSize::Get(attrs));

    return fit_map;
};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_lexicase = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return LexicaseFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_tournament = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return TournamentFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_eco_ea = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return EcoEAFitness(pop,attrs);};

template <typename PHEN_T>
std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> do_sharing = [](emp::vector<PHEN_T> & pop, all_attrs attrs=DEFAULT){return SharingFitness(pop,attrs);};


template <typename PHEN_T>
emp::WeightedGraph CalcCompetition(emp::vector<PHEN_T> pop, 
                        std::function<emp::vector<double>(emp::vector<PHEN_T>&, all_attrs)> fit_fun,
                        all_attrs attrs=DEFAULT) {

    emp::WeightedGraph effects(pop.size());


    emp::vector<double> fitnesses = fit_fun(pop, attrs);

    for (size_t i = 0; i < pop.size(); i++) {
        effects.SetLabel(i, emp::to_string(pop[i]));
        // std::cout << effects.GetLabel(i) << std::endl;

        emp::vector<PHEN_T> curr = pop;
        for (size_t ax = 0; ax < curr[i].size(); ax++) {
            curr[i][ax] = 0; // Replace org with null org so pop size stays same
        }

        emp::vector<double> new_fits = fit_fun(curr, attrs);
        for (size_t j = 0; j < pop.size(); j++ ) {
            if (i == j) {continue;}

            // In terms of floating-point imprecision issues, it's much better
            // to check for equality than doing the subtraction and checking
            // for the result to be 0
            if (!almost_equal(fitnesses[j], new_fits[j], 10)) {
                effects.AddEdge(i, j, fitnesses[j] - new_fits[j]);
            }
            // std::cout << effect << std::endl;
        }
    }

    return effects;
}


#endif
