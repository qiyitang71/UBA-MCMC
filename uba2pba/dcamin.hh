#pragma once

#include <spot/parseaut/public.hh>
#include <spot/twaalgos/isunamb.hh>
#include <spot/twaalgos/hoa.hh>
#include <spot/twaalgos/sccfilter.hh>
#include <spot/twaalgos/determinize.hh>
#include <spot/twaalgos/translate.hh>
#include <spot/twaalgos/isdet.hh>
#include <spot/twaalgos/dualize.hh>
#include <spot/twa/bddprint.hh>
#include <spot/misc/version.hh>
#include <spot/misc/optionmap.hh>
#include <spot/twaalgos/split.hh>
#include <spot/twaalgos/complement.hh>
#include <spot/twaalgos/simulation.hh>
#include <spot/tl/parse.hh>
#include <spot/tl/print.hh>
#include <spot/tl/simplify.hh>

#include <unordered_map>
// #include <unordered_set>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>
#include <stack>
#include <set>
#include <list>
#include <cassert>
#include <utility>
// #include <boost/functional/hash.hpp> 


bool is_included(const spot::twa_graph_ptr &aut_a, const spot::twa_graph_ptr &aut_b);

bool is_equivalent(const spot::twa_graph_ptr &aut_a, const spot::twa_graph_ptr &aut_b);

void print_set(const std::set<unsigned> &S);

typedef std::pair<unsigned, unsigned> state_pair;

struct pair_hash
{
    std::size_t operator()(const state_pair &pair) const noexcept {
        return ((std::size_t)pair.first) << 32 | ((std::size_t)pair.second);
    }
};

// struct pair_equal
// {
//     bool operator()(const state_pair& p1, const state_pair& p2) const {
//         return p1.first == p2.first && p1.second == p2.second;
//     }
// };

class dca_minimiser
{
private:
    // ======= private data
    const spot::twa_graph_ptr dca;

    // std::unordered_map<state_pair, bool, boost::hash<state_pair>> langeq_relation;
    // std::unordered_map<state_pair, bool, boost::hash<state_pair>> subsafe_relation;
    std::unordered_map<state_pair, bool, pair_hash> langeq_relation;
    std::unordered_map<state_pair, bool, pair_hash> subsafe_relation;

    std::vector<std::set<unsigned>> safe_components;
    std::vector<unsigned> safe_scc_info;
    std::set<unsigned> centralised_states;

    spot::twa_graph_ptr norm_aut_a;
    spot::twa_graph_ptr norm_aut_b;

    spot::twa_graph_ptr norm_rev_aut;

    spot::twa_graph_ptr safe_aut_a;
    spot::twa_graph_ptr safe_aut_b;

    // construct (safe) automaton for checking language inclusion
    // spot::twa_graph_ptr
    // get_aut(const spot::twa_graph_ptr aut, bool keep_acc_edge);

    spot::twa_graph_ptr
    get_aut(const spot::twa_graph_ptr aut, bool keep_acc_edge, bool rev_trans = false);

    state_pair 
    make_ord_pair(unsigned p, unsigned q)
    {
        if (p > q) {
            std::swap(p, q);
        }
        // smaller one first
        return std::make_pair(p, q);
    }
    bool
    language_equivalent(unsigned p, unsigned q)
    {
        if (p == q)
            return true;
        // smaller one is first
        state_pair key = make_ord_pair(p, q);
        auto it = langeq_relation.find(key);

        if (it != langeq_relation.end()) {
            return it->second;
        }else {
            bool is_eq = check_equivalence_opt(p, q);
            // std::cout << "langeq " << p << " " << q << " " << is_eq << std::endl;
            langeq_relation.emplace(key, is_eq);
            return is_eq;
        }
    }

    bool check_equivalence(unsigned p, unsigned q) 
    {
        // std::cout << "Check equivalence " << norm_aut_a << ", " << norm_aut_b << std::endl;
        norm_aut_a->set_init_state(p);
        norm_aut_b->set_init_state(q);
        return is_equivalent(norm_aut_a, norm_aut_b);
    }

    bool check_equivalence_opt(unsigned p, unsigned q); 

    void propagate_eq_result(const spot::twa_graph_ptr aut_eq, unsigned p, unsigned q, bool eq_result);

    bool
    strongly_equivalent(unsigned p, unsigned q)
    {
        if (p == q) return true;
        state_pair ord_key = make_ord_pair(p, q);
        auto it = langeq_relation.find(ord_key); 
        if ( it != langeq_relation.end() && ! it->second) {
            return false;
        }
        state_pair key = std::make_pair(p, q);
        it = subsafe_relation.find(key);
        if (it != subsafe_relation.end() && ! it->second) {
            return false;
        }
        key = std::make_pair(q, p);
        it = subsafe_relation.find(key);
        if (it != subsafe_relation.end() && ! it->second) {
            return false;
        }
        return subsafe_contain(p, q) && subsafe_contain(q, p) && language_equivalent(p, q);
    }

    bool check_subsafe(unsigned p, unsigned q) 
    {
        safe_aut_a->set_init_state(p);
        safe_aut_b->set_init_state(q);
        return is_included(safe_aut_a, safe_aut_b);
    }

    bool check_subsafe2(unsigned p, unsigned q);
    // void backtrack(const spot::twa_graph_ptr aut, int state);

    bool
    subsafe_contain(unsigned p, unsigned q)
    {
        if (p == q) return true;
        state_pair key = std::make_pair(p, q);
        auto it = subsafe_relation.find(key);
        if (it != subsafe_relation.end()) {
            return it->second;
        }else {
            bool is_subsafe = check_subsafe2(p, q);
            // bool is_subsafe = check_subsafe(p, q);
            // subsafe_relation.emplace(key, is_subsafe);
            // std::cout << "added " << p << " " << q << " "<< is_subsafe << std::endl;
            return is_subsafe;
        }
    }

    bool
    subsafe_equivalent(unsigned p, unsigned q)
    {
        if (p == q) return true;
        state_pair key = std::make_pair(p, q);
        state_pair ord_key = make_ord_pair(p, q);
        auto it = langeq_relation.find(ord_key); 
        if ( it != langeq_relation.end() && ! it->second) {
            return false;
        }
        it = subsafe_relation.find(key);
        if (it != subsafe_relation.end() && ! it->second) {
            return false;
        }

        return language_equivalent(p, q) && subsafe_contain(p, q);
    }

    // we need to pick the set of components that H(S, S')
    // p \in S, q \in S' such that L(p) = L(q) & Ls(p) <= Ls(q)
    bool
    safe_H_relation(const std::set<unsigned> &S, const std::set<unsigned> &Sp)
    {
        auto p = (*S.begin());
        for (const auto &q : Sp)
        {
            if (subsafe_equivalent(p, q))
            {
                return true;
            }
            // std::cout << false << std::endl;
        }
        return false;
    }

    // Normalization step: adjust transitions between components
    // that is, transitions between different safe components must be rejecting
    spot::twa_graph_ptr
    normalize(const spot::twa_graph_ptr &aut);

    // Safe centralization
    // q subsafe s, then they both belong to the same SCC
    spot::twa_graph_ptr
    safe_centralize(const spot::twa_graph_ptr &aut);

    // Minimize automaton using equivalence classes
    // q strongly equivalent to s, then they should be merged together
    spot::twa_graph_ptr
    safe_minimize(const spot::twa_graph_ptr &aut);

    void 
    prepare_subsafe_relation(const spot::twa_graph_ptr& aut);

    void
    compute_safe_components();

    void 
    prepare_langeq_relation(const spot::twa_graph_ptr& aut); 

public:

    // Precondition: aut does not have dead states
    //, which means all states are reachable
    dca_minimiser(const spot::twa_graph_ptr &aut) : dca(aut), safe_scc_info(aut->num_states(), -1)
    {
        if (!spot::is_deterministic(aut) || aut->acc().get_acceptance() != spot::acc_cond::acc_code::cobuchi()) {
            std::runtime_error("Automaton is not a DCA");
        }
        std::cout << "Computing safe components" << std::endl;
        compute_safe_components(); 
        std::cout << "Computing language equivalence" << std::endl;
        prepare_langeq_relation(dca);
    }

    // Main minimization function
    spot::twa_graph_ptr minimize()
    {
        std::cout << "Now normalise DCA" << std::endl;
        spot::twa_graph_ptr aut = normalize(dca);
        // spot::print_hoa(std::cout, aut);
        // std::cout << std::endl;
        // we need to make sure that aut is normal before
        // computing the subsafe relation
        std::cout << "Computing subsafe language equivalence" << std::endl;
        prepare_subsafe_relation(aut);
        std::cout << "Now safe centralise DCA" << std::endl;
        aut = safe_centralize(aut);
        std::cout << "Now safe minimise DCA" << std::endl;
        aut = safe_minimize(aut);
        return aut;
    }

};