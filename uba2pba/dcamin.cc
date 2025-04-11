#include "dcamin.hh"

using namespace spot;

bool is_included(const spot::twa_graph_ptr &aut_a, const spot::twa_graph_ptr &aut_b)
{
    auto complement_aut_b = spot::complement(aut_b);
    return !aut_a->intersects(complement_aut_b);
}

bool is_equivalent(const spot::twa_graph_ptr &aut_a, const spot::twa_graph_ptr &aut_b)
{
    return is_included(aut_a, aut_b) && is_included(aut_b, aut_a);
}

void print_set(const std::set<unsigned> &S)
{
    std::cout << "{";
    for (const auto &s : S)
    {
        std::cout << " " << s;
    }
    std::cout << "}";
}

// ======= private functions
// spot::twa_graph_ptr
// dca_minimiser::get_aut(const spot::twa_graph_ptr aut, bool keep_acc_edge)
// {
//     // We need to use the same dictionary
//     // and the same state numbers
//     // spot::twa_graph_ptr result = make_twa_graph(aut->get_dict());
//     // result->copy_ap_of(aut);
//     // result->new_states(aut->num_states());
//     // result->set_acceptance(1, spot::acc_cond::fin({0}));
//     // // result->set_init_state(init_state);
//     // for (unsigned src = 0; src < aut->num_states(); src++)
//     // {
//     //     for (auto &edge : aut->out(src))
//     //     {
//     //         bool acc = edge.acc.has(0);
//     //         if (!acc)
//     //         {
//     //             result->new_edge(src, edge.dst, edge.cond);
//     //         }
//     //         else if (keep_acc_edge && acc)
//     //         {
//     //             result->new_acc_edge(src, edge.dst, edge.cond);
//     //         }
//     //     }
//     // }
//     // std::cout << "automaton " << init_state << " acc = " << keep_acc_edge << std::endl;
//     // spot::print_hoa(std::cout, result);
//     return get_aut(aut, keep_acc_edge, false);
// }

spot::twa_graph_ptr
dca_minimiser::get_aut(const spot::twa_graph_ptr aut, bool keep_acc_edge, bool rev_trans)
{
    // We need to use the same dictionary
    // and the same state numbers
    spot::twa_graph_ptr result = make_twa_graph(aut->get_dict());
    result->copy_ap_of(aut);
    result->new_states(aut->num_states());
    result->set_acceptance(1, spot::acc_cond::fin({0}));
    // result->set_init_state(init_state);
    for (unsigned src = 0; src < aut->num_states(); src++)
    {
        for (auto &edge : aut->out(src))
        {
            bool acc = edge.acc.has(0);
            unsigned left_state = rev_trans ? edge.dst : src;
            unsigned right_state = rev_trans ? src : edge.dst;
            if (!acc)
            {
                result->new_edge(left_state, right_state, edge.cond);
            }
            else if (keep_acc_edge && acc)
            {
                result->new_acc_edge(left_state, right_state, edge.cond);
            }
        }
    }
    // std::cout << "automaton " << init_state << " acc = " << keep_acc_edge << std::endl;
    // spot::print_hoa(std::cout, result);
    return result;
}
// Normalization step: adjust transitions between components
// that is, transitions between different safe components must be rejecting
twa_graph_ptr
dca_minimiser::normalize(const twa_graph_ptr &aut)
{
    spot::twa_graph_ptr result = make_twa_graph(aut->get_dict());
    result->copy_ap_of(aut);
    result->new_states(aut->num_states());
    result->set_acceptance(1, spot::acc_cond::fin({0}));
    result->set_init_state(aut->get_init_state_number());
    for (auto src = 0; src < aut->num_states(); ++src)
    {
        auto src_scc = safe_scc_info[src];
        for (const auto &edge : aut->out(src))
        {
            auto dst_scc = safe_scc_info[edge.dst];
            if (dst_scc != src_scc || edge.acc.has(0))
            {
                // Adjust transition to rejecting if crossing components
                result->new_acc_edge(src, edge.dst, edge.cond);
            }
            else
            {
                result->new_edge(src, edge.dst, edge.cond);
            }
        }
    }
    // assert (is_equivalent(result, aut)) : "Normlise not right" ;
    return result;
}

// Safe centralization
// q subsafe s, then they both belong to the same SCC
twa_graph_ptr
dca_minimiser::safe_centralize(const twa_graph_ptr &aut)
{
    // we need to pick the set of components that H(S, S')
    // p \in S, q \in S' such that L(p) = L(q) & Ls(p) <= Ls(q)

    // We now compute the frontier F \subset SafeComponents that
    // for every safe component S, there exists S' \in F such that H(S, S')
    // and for two S, S' \in F, ¬H(S, S') and ¬H(S', S)

    std::list<std::set<unsigned>> frontier;
    for (auto &S : safe_components)
    {
        // std::cout << "safe component S: ";
        // print_set(S);
        // std::cout << std::endl;
        if (frontier.size() <= 0)
        {
            frontier.emplace_back(S);
            continue;
        }
        // decide whether it can be added or not
        std::list<std::set<unsigned>>::iterator it;
        bool replace_old = false;
        for (it = frontier.begin(); it != frontier.end(); ++it)
        {
            // std::cout << "safe component Sp: ";
            // print_set(*it);
            // std::cout << std::endl;
            bool cover_S = safe_H_relation(S, *it);
            if (cover_S)
            {
                // S is covered by current element
                break;
            }
            // This implies not H(S, S')
            // if H(S', S) and not H(S, S'), then remove S'
            bool cover_it = safe_H_relation(*it, S);
            if (cover_it)
            {
                // std::cout << "Covered" << std::endl;
                replace_old = true;
                break;
            }
            // else S' is not covered by S
        }

        if (replace_old)
        {
            frontier.erase(it);
            frontier.emplace_back(S);
            // std::cout << "Added ";
            // print_set(S);
            // std::cout << std::endl;
        }
        else if (it == frontier.end())
        {
            // by default, add this
            frontier.emplace_back(S);
            // std::cout << "Added ";
            // print_set(S);
            // std::cout << std::endl;
        }
    }
    // are there safe components S, S' such that H(S, S') and H(S', S)?
    // we need to careful
    // std::cout << "Computed frontier " << std::endl;
    for (const auto &safe_scc : frontier)
    {
        // std::cout << "safe scc: " ;
        // print_set(safe_scc);
        // std::cout << std::endl;
        centralised_states.insert(safe_scc.begin(), safe_scc.end());
    }
    // std::cout << "Remaining states:";
    // print_set(centralised_states);
    // std::cout << std::endl;
    std::vector<bool> remaining_states(aut->num_states(), false);
    for (const auto &s : centralised_states)
    {
        remaining_states[s] = true;
    }

    // Find initial state in safe components
    int initial_state = -1;
    unsigned old_init_number = aut->get_init_state_number();

    for (const auto &S : frontier)
    {
        if (S.find(old_init_number) != S.end())
        {
            initial_state = old_init_number;
            break;
        }
    }
    // std::cout << "Init: " << old_init_number << std::endl;

    if (initial_state == -1)
    {
        for (const auto &S : frontier)
        {
            bool out = false;
            for (const auto &state : S)
            {
                if (subsafe_equivalent(old_init_number, state))
                {
                    initial_state = state;
                    out = true;
                    break;
                }
            }
            if (out)
                break;
        }
    }
    // std::cout << "Init: " << initial_state << std::endl;

    // now we construct new automaton
    spot::twa_graph_ptr result = make_twa_graph(aut->get_dict());
    result->copy_ap_of(aut);
    result->new_states(aut->num_states());
    result->set_acceptance(1, spot::acc_cond::fin({0}));
    result->set_init_state(initial_state);
    for (auto src = 0; src < aut->num_states(); ++src)
    {
        if (!remaining_states[src])
        {
            continue;
        }
        // std::set<std::pair<bdd, unsigned>> acc_edges;
        // std::set<unsigned> rej_succs;
        // our automaton is still deterministic
        for (const auto &edge : aut->out(src))
        {
            // since it is still deterministic
            if (edge.acc.has(0))
            {
                // we redirect this transition to a state in different safe component
                // successors are those states equivalent to edge.dst
                // in remaining_states
                std::set<unsigned> succs;
                for (const auto &s : centralised_states)
                {
                    if (language_equivalent(s, edge.dst))
                    {
                        result->new_acc_edge(src, s, edge.cond);
                    }
                }
            }
            else
            {
                // if safe transition, then keep it
                // for sure edge.dst is in the same safe component
                // since input aut is normal
                result->new_edge(src, edge.dst, edge.cond);
            }
        }
    }
    // std::cout << "Centralised Equiv: " << is_equivalent(result, aut) << std::endl;
    return result;
}

// Minimize automaton using equivalence classes
// q strongly equivalent to s, then they should be merged together
twa_graph_ptr
dca_minimiser::safe_minimize(const twa_graph_ptr &aut)
{
    // we only consider merge equivalent states in remaining_states
    // std::vector<std::set<unsigned>> equivalence_classes;
    // compute the equivalence classes according to strongly equivalence
    std::map<unsigned, unsigned> eq_map;
    for (const auto &s : centralised_states)
    {
        bool found = false;
        for (unsigned t = 0; t < s; t++)
        {
            if (centralised_states.find(t) != centralised_states.end() && strongly_equivalent(s, t))
            {
                eq_map[s] = t;
                found = true;
                break;
            }
        }
        if (!found)
        {
            eq_map[s] = s;
        }
    }
    // equivalence classes print
    // for (const auto &item : eq_map)
    // {
    //     std::cout << item.first << " -> " << item.second << std::endl;
    // }
    // now we construct new automaton
    spot::twa_graph_ptr result = make_twa_graph(aut->get_dict());
    result->copy_ap_of(aut);
    result->new_states(aut->num_states());
    result->set_acceptance(1, spot::acc_cond::fin({0}));
    result->set_init_state(eq_map[aut->get_init_state_number()]);
    for (const auto &src : centralised_states)
    {
        for (const auto &edge : aut->out(src))
        {
            if (edge.acc.has(0))
            {
                result->new_acc_edge(eq_map[src], eq_map[edge.dst], edge.cond);
            }
            else
            {
                result->new_edge(eq_map[src], eq_map[edge.dst], edge.cond);
            }
        }
    }
    // std::cout << "Min Equiv: " << is_equivalent(result, aut) << std::endl;
    result->purge_dead_states();
    return result;
}

void dca_minimiser::prepare_subsafe_relation(const twa_graph_ptr &aut)
{
    safe_aut_a = get_aut(aut, false);
    safe_aut_b = get_aut(aut, false);
    // we need to store both direction
    // subsafe_relation.clear();
    // for (unsigned s = 0; s < aut->num_states(); s++)
    // {
    //     std::vector<bool> relation(aut->num_states(), false);
    //     safe_aut_a->set_init_state(s);
    //     for (unsigned t = 0; t < aut->num_states(); t++)
    //     {
    //         if (s == t)
    //         {
    //             relation[t] = true;
    //             continue;
    //         }
    //         safe_aut_b->set_init_state(t);
    //         bool subsafe = is_included(safe_aut_a, safe_aut_b);
    //         // std::cout << "SafeLang inclusion: L(" << s << ")  < L(" << t << "): " << subsafe << std::endl;
    //         // bool langequiv = s < t ? langequiv_relation[t][s] : langequiv_relation[s][t];
    //         relation[t] = subsafe;
    //     }
    //     subsafe_relation.emplace_back(relation);
    // }
    // for (unsigned s = 0; s < aut->num_states(); s++)
    // {
    //     for (unsigned t = 0; t < aut->num_states(); t++)
    //     {
    //         std::cout << "SafeLang : L(" << s << ")  < L(" << t << "): " << subsafe_relation[s][t] << std::endl;
    //     }
    // }
}

void dca_minimiser::compute_safe_components()
{
    auto safe_aut_a = get_aut(dca, false);
    std::vector<bool> visited(safe_aut_a->num_states(), false);
    // find all safe components
    for (unsigned s = 0; s < safe_aut_a->num_states(); s++)
    {
        if (visited[s])
        {
            continue;
        }
        // std::cout << "Init: " << s << std::endl;
        safe_aut_a->set_init_state(s);
        // now compute the safe components
        auto si = spot::scc_info(safe_aut_a);
        // std::cout << "scc num: " << si.scc_count() << std::endl;
        for (unsigned scc = 0; scc < si.scc_count(); scc++)
        {
            auto &states = si.states_of(scc);
            for (auto &state : states)
            {
                // std::cout << " visited " << state << std::endl;
                visited[state] = true;
            }
            std::set<unsigned> scc_states(states.begin(), states.end());
            auto it = std::find(safe_components.begin(), safe_components.end(), scc_states);
            if (it == safe_components.end())
                safe_components.emplace_back(scc_states);
        }
    }

    for (unsigned scc = 0; scc < safe_components.size(); scc++)
    {
        // std::cout << "safe scc " << scc << " : ";
        for (const auto &state : safe_components[scc])
        {
            safe_scc_info[state] = scc;
            // std::cout << " " << state;
        }
        // std::cout << std::endl;
    }
}

void dca_minimiser::prepare_langeq_relation(const twa_graph_ptr &aut)
{
    norm_aut_a = get_aut(aut, true);
    norm_aut_b = get_aut(aut, true);
    norm_rev_aut = get_aut(aut, true, true);
    // langeq_relation.clear();
    // for (unsigned s = 0; s < aut->num_states(); s++)
    // {
    //     std::vector<bool> relation(s, false);
    //     norm_aut_a->set_init_state(s);
    //     for (unsigned t = 0; t < s; t++)
    //     {
    //         norm_aut_b->set_init_state(t);
    //         bool lang_equiv = is_included(norm_aut_a, norm_aut_b) && is_included(norm_aut_b, norm_aut_a);
    //         relation[t] = lang_equiv;
    //         // std::cout << "Lang equiv: L(" << s << ")  = L(" << t << "): " << lang_equiv << std::endl;
    //     }
    //     std::cout << "Finished state " << s << std::endl;
    //     langeq_relation.emplace_back(relation);
    // }
}
/// @brief whether the safe language of p is included in safe language of q
/// we compute a fixed set of state pairs, if no state pair is disproved
/// then we already prove that they are safe
/// @param p
/// @param q
/// @return true if included
bool dca_minimiser::check_subsafe2(unsigned p, unsigned q)
{
    // std::cout << "check pair (" << p << ", " << q << ")" << std::endl;
    std::set<state_pair> visited;
    // visited.insert(std::make_pair(p, q));
    std::list<state_pair> todos;
    std::unordered_map<unsigned, std::set<unsigned>> parent_relation;
    // auto dict = spot::make_bdd_dict();
    // // check safe automata
    const auto p_scc_size = safe_components[safe_scc_info[p]].size();
    const auto q_scc_size = safe_components[safe_scc_info[q]].size();

    // product->set_acceptance(1, spot::acc_cond::fin({0}));
    // every state (x,y) in the product will be encoded as
    // q_scc_size * x + y
    auto get_state = [&q_scc_size](unsigned a, unsigned b)
    { return a * q_scc_size + b; };
    auto get_pair = [&q_scc_size](unsigned state)
    {
        unsigned b = state % q_scc_size;
        unsigned a = state / q_scc_size;
        return std::make_pair(a, b);
    };
    todos.push_back(std::make_pair(p, q));
    // unsigned num_added_parent = 0;
    while (todos.size() > 0)
    {
        // compute a fixed set of state pairs
        // should not be reference to front element
        auto curr = todos.front();
        todos.pop_front();

        // all accepting transitions have been removed
        // note that our automata are deterministic
        // we only check inclusion, not equivalence
        if (visited.find(curr) != visited.end())
        {
            continue;
        }
        // if (visited.find(curr) != visited.end()) {
        // keep adding their successor pairs
        unsigned curr_state = get_state(curr.first, curr.second);
        for (const auto &e1 : safe_aut_a->out(curr.first))
        {
            auto all_letters = e1.cond;
            for (const auto &e2 : safe_aut_b->out(curr.second))
            {
                auto succ = std::make_pair(e1.dst, e2.dst);
                auto letter = e1.cond & e2.cond;
                // std::cout << " succ pair: " << succ.first << ", " << succ.second << std::endl;
                // std::cout << " letter "<< letter << std::endl;
                // need to check (p', q') where p' != q'
                if (letter != bddfalse && e1.dst != e2.dst )
                {
                    // need to recursively check succ 
                    if (visited.find(succ) == visited.end()) {
                        todos.push_back(succ);
                    }
                    // std::cout << " succ pair: " << succ.first << ", " << succ.second << std::endl;
                    // std::cout << " curr pair: " << curr.first << ", " << curr.second << std::endl;
                    // successor point to parent
                    unsigned succ_state = get_state(e1.dst, e2.dst);

                    auto it = parent_relation.find(succ_state);
                    if (it == parent_relation.end()) {
                        std::set<unsigned> parents;
                        parents.insert(curr_state);
                        parent_relation.emplace(succ_state, parents);
                    }else {
                        it->second.insert(curr_state);
                    }
                }
                all_letters -= letter;
            }
            // if some letter is not covered
            if (all_letters != bddfalse)
            {
                // std::cout << " left edge: ";
                // spot::bdd_print_formula(std::cout, safe_aut_a->get_dict(), all_letters);
                // std::cout << std::endl;

                // not included, then all pairs were not included
                // language inclusion of current pair not satisfied
                // backtrack to all parents
                // auto state_id = get_state(curr.first, curr.second);
                std::set<unsigned> reach_states;
                // traverse all successors
                std::list<unsigned> work_list;
                work_list.push_back(curr_state);
                while (work_list.size() > 0 )
                {
                    // traverse all reachable states
                    unsigned state_id = work_list.front();
                    work_list.pop_front();
                    if (reach_states.find(state_id) != reach_states.end())
                    {
                        continue;
                    }
                    reach_states.insert(state_id);
                    auto it = parent_relation.find(state_id);
                    if (it == parent_relation.end()) {
                        // we already reach the initial pair
                        continue;
                    }
                    for (const auto &s : it->second)
                    {
                        if (reach_states.find(s) != reach_states.end())
                        {
                            continue;
                        }
                        work_list.push_back(s);
                    }
                    // break;
                }
                
                // now we set all states
                for (const auto state : reach_states)
                {
                    auto pair = get_pair(state);
                    subsafe_relation.emplace(pair, false);
                }
                return false;
            }
        }
        visited.insert(curr);
    }
    // all pairs should be language included
    for (const auto& pair : visited)
    {
        // std::cout << "true pair: " << pair.first << ", " << pair.second << std::endl;
        subsafe_relation.emplace(pair, true);
        // subsafe_relation[new_pair] =  true;
    }
    // std::cout << " return " << std::endl;
    return true;
}

/// @brief whether the language of p is equivalent to the language of q
/// we first use check equivalence 
/// @param p
/// @param q
/// @return true if equivalent
bool dca_minimiser::check_equivalence_opt(unsigned p, unsigned q)
{
    bool is_eq = check_equivalence(p, q);
    // we propagate the results 
    if (is_eq) {
        // propagate the successors
        propagate_eq_result(norm_aut_a, p, q, is_eq);
    }else {
        // propagate the predecessors
        propagate_eq_result(norm_rev_aut, p, q, is_eq);
    }
    return is_eq;
}

void dca_minimiser::propagate_eq_result(
    const spot::twa_graph_ptr aut_eq, unsigned p, unsigned q, bool eq_result)
{
        // std::cout << "check pair (" << p << ", " << q << ")" << std::endl;
    std::set<state_pair> visited;

    // visited.insert(std::make_pair(p, q));
    std::list<state_pair> todos;
    todos.push_back(make_ord_pair(p, q));
    // unsigned num_added_parent = 0;
    while (todos.size() > 0)
    {
        // compute a fixed set of state pairs
        // should not be reference to front element
        auto curr = todos.front();
        todos.pop_front();

        if (visited.find(curr) != visited.end())
        {
            continue;
        }
        // keep adding their successor pairs
        for (const auto &e1 : aut_eq->out(curr.first))
        {
            for (const auto &e2 : aut_eq->out(curr.second))
            {
                auto succ = make_ord_pair(e1.dst, e2.dst);
                auto letter = e1.cond & e2.cond;
                // std::cout << " succ pair: " << succ.first << ", " << succ.second << std::endl;
                // std::cout << " letter "<< letter << std::endl;
                // need to check (p', q') where p' != q'
                if (letter != bddfalse && e1.dst != e2.dst )
                {
                    // need to recursively check succ 
                    if (visited.find(succ) == visited.end()) {
                        todos.push_back(succ);
                    }
                    // std::cout << " succ pair: " << succ.first << ", " << succ.second << std::endl;
                    // std::cout << " curr pair: " << curr.first << ", " << curr.second << std::endl;
                    // successor point to parent
                    langeq_relation.emplace(succ, eq_result);
                }
            }
        }
        visited.insert(curr);
    }

}


/*
// Example usage
int main(int argc, char *argv[])
{

    std::string file_name;
    int verbose = 0;
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];

        if (arg == "-a")
        {
            if (argc < i + 1)
            {
                std::cerr << "ltl2pba: Option -a requires an argument.\n";
                return 1;
            }
            else
            {
                std::string str(argv[i + 1]);
                file_name = str;
                i++;
            }
        }else if (arg == "-v") {
            if (argc < i + 1)
            {
                std::cerr << "ltl2pba: Option -v requires an argument.\n";
                return 1;
            }
            else
            {
                std::string str(argv[i + 1]);
                verbose = std::stoi(str);
                i++;
            }
        }
    }

    auto dict = spot::make_bdd_dict();
    spot::automaton_stream_parser parser(file_name);
    spot::parsed_aut_ptr parsed_aut = parser.parse(dict);
    if (parsed_aut->format_errors(std::cerr))
    {
        std::runtime_error("File " + file_name + " is not in valid HOA format");
        return 1;
    }
    auto aut = parsed_aut->aut;

    if (!spot::is_deterministic(aut) || aut->get_acceptance() != spot::acc_cond::acc_code::cobuchi())
    {
        throw std::runtime_error("Accept only DCA");
    }

    // for (unsigned s = 0; s < aut->num_states(); s++)
    // {
    //     for (const auto &edge : aut->out(s))
    //     {
    //         std::cout << edge.cond << " " << edge.dst << std::endl;
    //     }
    // }

    aut->purge_dead_states();
    dca_minimiser dcamin(aut);
    auto res = dcamin.minimize();
    res = spot::split_edges(res);
    spot::print_hoa(std::cout, res);
    std::cout << std::endl;
    std::cout << "#states: " << res->num_states() << std::endl;
    std::cout << "Included: " << is_included(res, aut) << std::endl;
    std::cout << "Included: " << is_included(aut, res) << std::endl;
    return 0;
}
*/