/*********************************************************************
 *
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2008, 2013, Willow Garage, Inc.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of Willow Garage, Inc. nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 * Author: Eitan Marder-Eppstein
 *         David V. Lu!!
 *********************************************************************/
#include<global_planner/astar.h>
#include<costmap_2d/cost_values.h>

namespace global_planner {

AStarExpansion::AStarExpansion(PotentialCalculator* p_calc, int xs, int ys) :
        Expander(p_calc, xs, ys) {
}

bool AStarExpansion::calculatePotentials(unsigned char* costs, double start_x, double start_y, double end_x, double end_y,
                                        int cycles, float* potential) {
    //参数：costs地图指针、起点坐标、终点坐标、循环次数（地图栅格的二倍）、代价数组

    //下面是A*算法过程，使用队列 queue ,按照一定的优先级，将起点放入队列；将potential数组全部设为最大POT_HIGH；将起点的potential设为0
    queue_.clear();
    int start_i = toIndex(start_x, start_y);    //起点的索引
    queue_.push_back(Index(start_i, 0));    //起点加入队列

    std::fill(potential, potential + ns_, POT_HIGH);    //将所有网格的势能值初始化为一个极大值
    potential[start_i] = 0;

    int goal_i = toIndex(end_x, end_y);     //目标点索引
    int cycle = 0;
    //进入循环，只要队列不为空，得到最小cost的索引，并删除
    while (queue_.size() > 0 && cycle < cycles) {
        Index top = queue_[0];  //得到最小cost的索引并删除
        std::pop_heap(queue_.begin(), queue_.end(), greater1());    
        queue_.pop_back();

        int i = top.i;
        if (i == goal_i)
            return true;
        //将代价最小点i周围点加入搜索队里并更新代价值, 即对前后左右四个点执行add函数
        add(costs, potential, potential[i], i + 1, end_x, end_y);
        add(costs, potential, potential[i], i - 1, end_x, end_y);
        add(costs, potential, potential[i], i + nx_, end_x, end_y);
        add(costs, potential, potential[i], i - nx_, end_x, end_y);

        cycle++;
    }

    return false;
}

void AStarExpansion::add(unsigned char* costs, float* potential, float prev_potential, int next_i, int end_x,
                         int end_y) {
    //避免重复扩展和非法网格
    if (next_i < 0 || next_i >= ns_)
        return;

    if (potential[next_i] < POT_HIGH)
        return;

    if(costs[next_i]>=lethal_cost_ && !(unknown_ && costs[next_i]==costmap_2d::NO_INFORMATION))
        return;

    //计算节点距离起点的代价
    potential[next_i] = p_calc_->calculatePotential(potential, costs[next_i] + neutral_cost_, next_i, prev_potential);
    int x = next_i % nx_, y = next_i / nx_;
    float distance = abs(end_x - x) + abs(end_y - y);   // Manhattan距离
    
    //f(n) = g(n) + h(n)
    queue_.push_back(Index(next_i, potential[next_i] + distance * neutral_cost_));
    std::push_heap(queue_.begin(), queue_.end(), greater1());
}

} //end namespace global_planner
