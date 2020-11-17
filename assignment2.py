# Them cac thu vien neu can
import math
import random
import time
import tool

class Heuristic:
    def __init__(self, warehouse, N, M, order_lst, local_trials, global_trials):
        self.warehouse = warehouse # (xcoord, ycoord) of warehouse (A)
        self.N = N  # number of order
        self.M = M  # number of shipper
        self.order_lst = order_lst  # list of order, each order is a list of information: [xcoord, ycoord, volumn, weight]
        self.local_trials = local_trials # maximum number of continous non-better next states in one local hill climbing 
        self.global_trials = global_trials # maximum number of continuous non-better local maxima
        self.tasks_table = self.init_tasks_table() # list of shipper's route, each route is a list of order id
        self.heuristic = self.compute_heuristic() # heuristic function -> the target is to minimize this value

    def init_tasks_table(self):
        orders = list(range(self.N))
        random.shuffle(orders)
        return [orders[i::self.M] for i in range(self.M)]

    def print_state(self):
        print(self.tasks_table, end=' ')
        print(self.heuristic)

    # Random restart Hill Climbing
    def solve(self):
        global_trials = self.global_trials
        current_best_state = self.tasks_table
        current_best_heuristic = self.heuristic
        while global_trials > 0:
            local_trials = self.local_trials
            while local_trials != 0:
                operation = bool(random.getrandbits(1)) # random choose operation: intra_swap or inter_insert
                if operation:
                    result = self.intra_swap()
                    if result:
                        local_trials = self.local_trials
                        self.print_state()
                    else:
                        local_trials -= 1
                else:
                    result = self.inter_insert()
                    if result:
                        local_trials = self.local_trials
                        self.print_state()
                    else:
                        local_trials -= 1    

            if current_best_heuristic > self.heuristic:
                current_best_heuristic = self.heuristic
                current_best_state = self.tasks_table
                global_trials = self.global_trials
            else:
                global_trials -= 1
            print("Finish") # finish one local search

            # random restart
            self.tasks_table = self.init_tasks_table()
            self.heuristic = self.compute_heuristic()

        self.tasks_table = current_best_state
        self.heuristic = current_best_heuristic

    # Swap 2 order within one shipper's route
    def intra_swap(self):
        target_shipper = random.randint(0, self.M-1)    # random pick a target shipper's route to perform swap
        # print(target_shipper)
        route = self.tasks_table[target_shipper] # old route
        prev_profit = self.f[target_shipper] # old route profit
        delta_heuristic = sum([abs(prev_profit - f_i) for f_i in self.f]) * 2
        if len(route) >= 2:
            swap_pos = random.sample(range(len(route)), k=2)
            # print(swap_pos)
            route[swap_pos[0]], route[swap_pos[1]] = route[swap_pos[1]], route[swap_pos[0]] # perform swap
            self.f[target_shipper] = self.profit(route) # update route profit

            delta_heuristic -= sum([abs(self.f[target_shipper] - f_i) for f_i in self.f]) * 2
            if delta_heuristic <= 0: 
                # Go back to previous state
                route[swap_pos[0]], route[swap_pos[1]] = route[swap_pos[1]], route[swap_pos[0]] # perform swap
                self.f[target_shipper] = self.profit(route) # update route profit
                assert self.f[target_shipper] == prev_profit
                return False                
            self.heuristic -= delta_heuristic
            return True
        return False

    # Move one order from one shipper (giver) to another shipper (receiver) 
    def inter_insert(self):
        giver, receiver = random.sample(range(self.M), k=2) # random pick giver and receiver
        giver_route = self.tasks_table[giver] # old route
        receiver_route = self.tasks_table[receiver] # old route
        prev_giver_profit = self.f[giver] # old route profit
        prev_receiver_profit = self.f[receiver] # old route profit
        delta_heuristic = ( sum([abs(prev_giver_profit - f_i) for f_i in self.f]) + sum([abs(prev_receiver_profit - f_i) for f_i in self.f]) - abs(prev_giver_profit - prev_receiver_profit) ) * 2
        if len(giver_route) >= 2:
            target_order_pos = random.randint(0, len(giver_route) - 1)
            target_order = giver_route[target_order_pos]
            insert_pos = random.randint(0, len(receiver_route))
            # print(swap_pos)
            # giver_route.remove(target_order) # remove target order from giver
            del giver_route[target_order_pos] # remove target order from giver
            receiver_route.insert(insert_pos, target_order) # insert target order to receiver at target position
            self.f[giver] = self.profit(giver_route) # update route profit
            self.f[receiver] = self.profit(receiver_route) # update route profit

            delta_heuristic -= ( sum([abs(self.f[giver] - f_i) for f_i in self.f]) + sum([abs(self.f[receiver] - f_i) for f_i in self.f]) - abs(self.f[giver] - self.f[receiver]) ) * 2
            if delta_heuristic <= 0: 
                # Go back to previous state
                # receiver_route.remove(target_order) # remove target order from receiver
                del receiver_route[insert_pos] # remove target order from receiver
                giver_route.insert(target_order_pos, target_order) # insert target order to giver at target position
                self.f[giver] = self.profit(giver_route) # update route profit
                self.f[receiver] = self.profit(receiver_route) # update route profit

                assert self.f[giver] == prev_giver_profit
                assert self.f[receiver] == prev_receiver_profit
                return False                
            self.heuristic -= delta_heuristic
            # print(self.heuristic)
            # print(self.compute_heuristic())
            # assert self.heuristic == self.compute_heuristic(), "abc"
            return True
        return False        

    # Calculate profit of a route
    def profit(self, route):
        income = 0
        distance = 0
        start = self.warehouse
        for order_id in route:
            order = self.order_lst[order_id]
            income += 5 + order[2] + 2*order[3]
            distance += math.sqrt( (start[0]-order[0])**2 + (start[1]-order[1])**2 )
            start = (order[0], order[1])
        return income - (distance / 40 * 20 + 10)

    # Target is to minimize this function
    def compute_heuristic(self):
        self.f = [self.profit(route) for route in self.tasks_table]
        return sum([abs(f_i - f_j) for f_i in self.f for f_j in self.f])
        


def read_input(file_input):
    with open(file_input, 'r') as f:
        data = f.readlines()
        warehouse = (int(data[0].split()[0]), int(data[0].split()[1]))
        N, M = int(data[1].split()[0]), int(data[1].split()[1])
        order_lst = [data[order_id+2].split() for order_id in range(N)]
        order_lst = [tuple([int(order_info) for order_info in order]) for order in order_lst]
        return warehouse, N, M, order_lst

def write_output(shipper_lst, file_output):
    with open(file_output, 'w') as f:
        shipper_lst = [[str(order) for order in shipper] for shipper in shipper_lst]
        f.writelines([' '.join(shipper) + '\n' for shipper in shipper_lst])

def assign(file_input, file_output):
    # generate test input
    # tool.generate_test(20, 5, file_input) # random generate input test

    # read input
    warehouse, N, M, order_lst = read_input(file_input)
    solver = Heuristic(warehouse, N, M, order_lst, local_trials=500, global_trials=10)
    
    # run algorithm
    solver.print_state()
    start = time.time()
    solver.solve()
    end = time.time()
    solver.print_state()
    print(sum(solver.f) / M)
    print(solver.heuristic / (M**2))
    print(solver.f)
    print("Time elapsed: " + str(end - start))
    # write output
    write_output(solver.tasks_table, file_output)

    return


assign('input1.txt', 'output1.txt')
