#include <iostream>
#include <vector>

#include "CGL/vector2D.h"

#include "mass.h"
#include "rope.h"
#include "spring.h"

namespace CGL {

    Rope::Rope(Vector2D start, Vector2D end, int num_nodes, float node_mass, float k, vector<int> pinned_nodes)
    {
        // TODO (Part 1): Create a rope starting at `start`, ending at `end`, and containing `num_nodes` nodes.
        // class member
        // vector<Mass *> masses;
        // vector<Spring *> springs;
        for (int i = 0; i < num_nodes; ++i)
        {
            Vector2D distancePerNode = (end - start) / (num_nodes - 1);
            masses.push_back(new Mass(
                start + i * distancePerNode, node_mass, false
            ));
        }
        
        for (int i = 0; i < num_nodes - 1; ++i)
        {
            springs.push_back(new Spring(
                masses[i], masses[i + 1], k       
            ));
        }

    //    Comment-in this part when you implement the constructor
       for (auto &i : pinned_nodes) {
           masses[i]->pinned = true;
       }
    }

    void Rope::simulateEuler(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 2): Use Hooke's law to calculate the force on a node
            auto& m1 = s->m1;  // & !!!
            auto& m2 = s->m2;
            auto ks = s->k;
            auto res_len = s->rest_length;
            m1->forces = m1->forces - ks * (m1->position - m2->position).unit() * 
                         ((m1->position - m2->position).norm() - res_len);
            m2->forces = m2->forces + ks * (m1->position - m2->position).unit() * 
                         ((m1->position - m2->position).norm() - res_len);
            // Vector2D ab = s->m2->position - s->m1->position;
            // Vector2D f = s->k * (ab / ab.norm()) * (ab.norm() - s->rest_length);
            // s->m1->forces += f;
            // s->m2->forces -= f;
        }

        for (auto &m : masses)
        {
            float kd = 0.01f;
            if (!m->pinned)
            {
                // TODO (Part 2): Add the force due to gravity, then compute the new velocity and position
                m->forces += m->mass * gravity;
                // TODO (Part 2): Add global damping
                // m->force -= kd * v;
                if (true)
                {
                    m->forces -= kd * m->velocity;
                }

                if (false)
                {
                    m->position += m->velocity * delta_t;
                    m->velocity += m->forces / m->mass * delta_t;
                }
                else
                {
                    m->velocity += m->forces / m->mass * delta_t;
                    m->position += m->velocity * delta_t;
                }
            }

            // Reset all forces on each mass
            m->forces = Vector2D(0, 0);
        }
    }

    void Rope::simulateVerlet(float delta_t, Vector2D gravity)
    {
        for (auto &s : springs)
        {
            // TODO (Part 3): Simulate one timestep of the rope using explicit Verlet （solving constraints)
            auto& m1 = s->m1;  // & !!!
            auto& m2 = s->m2;
            auto ks = s->k;
            auto res_len = s->rest_length;
            m1->forces = m1->forces - ks * (m1->position - m2->position).unit() * 
                         ((m1->position - m2->position).norm() - res_len);
            m2->forces = m2->forces + ks * (m1->position - m2->position).unit() * 
                         ((m1->position - m2->position).norm() - res_len);
        }

        for (auto &m : masses)
        {
            if (!m->pinned)
            {
                // TODO (Part 3.1): Set the new position of the rope mass
                m->forces += m->mass * gravity;
                Vector2D lastPosition = m->position;
                // m->position = m->position + (m->position - m->last_position) +
                //               m->forces / m->mass * delta_t * delta_t;
                // m->last_position = lastPosition;
                // TODO (Part 4): Add global Verlet damping
                float damping_factor = 0.00005f;
                m->position = m->position + damping_factor * (m->position - m->last_position) +
                              m->forces / m->mass * delta_t * delta_t;
                m->last_position = lastPosition;
            }
            m->forces =  Vector2D(0,0);
        }
    }
}
