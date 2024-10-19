#include "scan_matching_skeleton/correspond.h"
#include "cmath"

using namespace std;
#define MAX_DIST 1000000.0

const int UP_SMALL = 0;
const int UP_BIG = 1;
const int DOWN_SMALL = 2;
const int DOWN_BIG = 3;

void getCorrespondence(vector<Point> &old_points, vector<Point> &trans_points, vector<Point> &points,
                       vector<vector<int>> &jump_table, vector<Correspondence> &c, float prob)
{
    // Clear previous correspondences
    c.clear();
    int last_best = -1;
    const int trans_size = trans_points.size();
    const int old_size = old_points.size();

    // Loop through each transformed point to find its correspondence in the old points
    for (int ind_trans = 0; ind_trans < min(old_size, trans_size); ++ind_trans)
    {
        // Variables to store the best match and its distance
        int best = -1;
        int second_best = -1;
        double best_dist = MAX_DIST;
        double second_best_dist = MAX_DIST;

        // Reference to the current transformed point
        Point &p_trans = trans_points[ind_trans];

        // Calculate the approximate starting index in the old points based on angle
        double angle_diff = p_trans.theta - old_points[0].theta;
        int start_index = static_cast<int>((angle_diff) * (old_size / (2.0 * M_PI)));
        start_index = std::max(0, std::min(start_index, old_size - 1));

        // If last_best is valid, use it as the starting point
        int we_start_at = (last_best != -1) ? (last_best + 1) : start_index;

        // Initialize indices for searching up and down
        int up = we_start_at;
        int down = we_start_at - 1;
        double last_dist_up = MAX_DIST;
        double last_dist_down = MAX_DIST;
        bool up_stopped = false;
        bool down_stopped = false;

        // Search until both directions have been stopped
        while (!up_stopped || !down_stopped)
        {
            // Determine whether to search up or down based on the last distances
            bool now_up = !up_stopped && (last_dist_up < last_dist_down);

            if (now_up)
            {
                // If we've reached the end of the points, stop searching upwards
                if (up >= old_size)
                {
                    up_stopped = true;
                    continue;
                }

                // Calculate the distance between the current transformed point and the point in the "up" direction
                last_dist_up = p_trans.distToPoint2(&old_points[up]);

                // If this point is better than the best found so far, update best
                if (last_dist_up < best_dist)
                {
                    second_best = best;
                    second_best_dist = best_dist;
                    best = up;
                    best_dist = last_dist_up;
                }
                else if (last_dist_up < second_best_dist)
                {
                    second_best = up;
                    second_best_dist = last_dist_up;
                }

                // Check if we can use an early stopping condition (if moving away from the start index)
                if (up > start_index)
                {
                    double delta_phi = old_points[up].theta - p_trans.theta;
                    double min_dist_up = std::sin(delta_phi) * p_trans.r;

                    if (min_dist_up * min_dist_up > best_dist)
                    {
                        // Stop searching upwards if we can't find a better match
                        up_stopped = true;
                        continue;
                    }

                    // Use jump table optimization if possible
                    up = (old_points[up].r < p_trans.r) ? jump_table[up][UP_BIG] : jump_table[up][UP_SMALL];
                }
                else
                {
                    // If moving toward the start index, increment without using jump table
                    up++;
                }
            }
            else
            {
                // If we've reached the start of the points, stop searching downwards
                if (down < 0)
                {
                    down_stopped = true;
                    continue;
                }

                // Calculate the distance between the current transformed point and the point in the "down" direction
                last_dist_down = p_trans.distToPoint2(&old_points[down]);

                // If this point is better than the best found so far, update best
                if (last_dist_down < best_dist)
                {
                    second_best = best;
                    second_best_dist = best_dist;
                    best = down;
                    best_dist = last_dist_down;
                }
                else if (last_dist_down < second_best_dist)
                {
                    second_best = down;
                    second_best_dist = last_dist_down;
                }

                // Check if we can use an early stopping condition (if moving away from the start index)
                if (down < start_index)
                {
                    double delta_phi = old_points[down].theta - p_trans.theta;
                    double min_dist_down = std::sin(delta_phi) * p_trans.r;

                    if (min_dist_down * min_dist_down > best_dist)
                    {
                        // Stop searching downwards if we can't find a better match
                        down_stopped = true;
                        continue;
                    }

                    // Use jump table optimization if possible
                    down = (old_points[down].r < p_trans.r) ? jump_table[down][DOWN_BIG] : jump_table[down][DOWN_SMALL];
                }
                else
                {
                    // If moving toward the start index, decrement without using jump table
                    down--;
                }
            }
        }

        // If a best match was found, create the correspondence, otherwise add a null correspondence
        if (best != -1 && second_best != -1)
        {
            c.push_back(Correspondence(&trans_points[ind_trans], &points[ind_trans], &old_points[best], &old_points[second_best]));
        }
        else
        {
            // If no valid correspondence is found, add a null correspondence
            c.push_back(Correspondence(&trans_points[ind_trans], &points[ind_trans], nullptr, nullptr));
        }

        // Update last_best to use as a hint in the next iteration
        last_best = best;
    }
}

void computeJump(vector<vector<int>> &table, vector<Point> &points)
{
  table.clear();
  int n = points.size();
  for (int i = 0; i < n; ++i)
  {
    vector<int> v = {n, n, -1, -1};
    for (int j = i + 1; j < n; ++j)
    {
      if (points[j].r < points[i].r)
      {
        v[UP_SMALL] = j;
        break;
      }
    }
    for (int j = i + 1; j < n; ++j)
    {
      if (points[j].r > points[i].r)
      {
        v[UP_BIG] = j;
        break;
      }
    }
    for (int j = i - 1; j >= 0; --j)
    {
      if (points[j].r < points[i].r)
      {
        v[DOWN_SMALL] = j;
        break;
      }
    }
    for (int j = i - 1; j >= 0; --j)
    {
      if (points[j].r > points[i].r)
      {
        v[DOWN_BIG] = j;
        break;
      }
    }
    table.push_back(v);
  }
}
