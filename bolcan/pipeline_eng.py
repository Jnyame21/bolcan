import math
import pandas as pd

def pipe_network(loops, equation:str, pipes_in_common_loops=None, K=False, error=1e-8, max_iter=100):
    """
    Find the rate and direction of flow in pipes in a pipe network.

    Parameters
    -----------
    loops: list, (required)
        A list of dictionaries that contains the name and pipes in that loop.
    
    equation: str, required
        The equation to use for the calculation (e.g., 'darcy' for Darcy-Weisbach or 'hazen' for Hazen-Williams).
    
    pipes_in_common_loops: list, default(None), (optional)
        A list containing dictionaries of pipes and the loops they are common to. Eg: [{'name': 'A-B', 'loops': ['Loop A', 'Loop B']}]
    
    K: boolean, default(False), (optional)
        If True, K values will be used instead of the friction factor (f) when using the Darcy-Weisbach equation or the Hazen-Williams coefficient (C) when using the Hazen-Williams equation.
    
    error: float, default(1e-8), (optional)
        The tolerance of the rate.
    
    max_iter: int, default(100), (optional)
        The maximum number of iterations.

    Returns
    --------
    A list of dictionaries containing each loop's name and a list of each pipe's data together with the assumed and corrected rates.
        
    Examples
    --------
    Without K(thus with K not specified or set to False):
        loops = [
            {'name': 'Loop A', 'props': [
                ['pipe', 'length', 'diameter', 'C', 'Qa'],
                ['A-B', 500, 0.2, 100, 1],
                ['A-C', 330, 0.35, 100, -11],
                ['B-D', 330, 0.2, 100, 0.5],
                ['C-D', 500, 0.2, 100, -1],
            ]}
        ]
        
        pipe_network(loops, equation='hazen')
        
        Result:
            [{'name': 'Loop A',
            'props': [
            ['pipe', 'length', 'diameter', 'C', 'Qa', 'Qnew', 'hl', 'hl/Qnew'],
            ['A-B', 500, 0.2, 100, 1.0, 1.6547, 6850.7357, 4140.2285],
            ['A-C', 330, 0.35, 100, -11.0, -10.3453, -8797.1821, 850.3534],
            ['B-D', 330, 0.2, 100, 0.5, 1.1547, 2323.8782, 2012.5808],
            ['C-D', 500, 0.2, 100, -1.0, -0.3453, -377.4318, 1092.9779]]
            }]
        
        
    With K(thus with K set to True):
        loops = [
            {'name': 'Loop A', 'props': [
                ['pipe', 'K', 'Qa'],
                ['A-B', 2, -20],
                ['A-D', 1, 80],
                ['B-D', 3, 20],
            ]},
            {'name': 'Loop B', 'props': [
                ['pipe', 'K', 'Qa'],
                ['B-D', 3, -20],
                ['B-C', 1, -40],
                ['C-D', 2, 35],
            ]},
        ]

        common_pipes = [
            {'name': 'B-D', 'loops': ['Loop A', 'Loop B']},
        ]
        
        pipe_network(loops, pipes_in_common_loops=common_pipes, equation='darcy', K=True)
        
        Result:
            [{'name': 'Loop A',
            'props': [
            ['pipe', 'K', 'Qa', 'Qnew', 'hl', 'hl/Qnew'],
            ['A-B', 2, -20, -41.4811, -3441.3607, 82.9622],
            ['A-D', 1, 80, 58.5189, 3424.4635, 58.5189],
            ['B-D', 3, 20, 2.3733, 16.8973, 7.1198]]},
            {'name': 'Loop B',
            'props': [
            ['pipe', 'K', 'Qa', 'Qnew', 'hl', 'hl/Qnew'],
            ['B-D', 3, -20, -2.3733, -16.8973, 7.1198],
            ['B-C', 1, -40, -43.8544, -1923.2047, 43.8544],
            ['C-D', 2, 35, 31.1456, 1940.102, 62.2913]]
            }]
        
    """

    if equation not in ['darcy', 'hazen']:
        raise ValueError("The equation must either be 'hazen' for Hazen Williams or 'darcy' for Darcy Weisbach")
    
    n = 2 if equation == 'darcy' else 1.85
    
    df_list = []
    for lp in loops:
        header = lp['props'][0]
        data = lp['props'][1:]
        df = pd.DataFrame(data, columns=header)
        
        if 'pipe' not in df.columns:
            raise ValueError("No 'pipe' column was found. Make sure you have the names of the pipes under a column named 'pipe' ")
        elif 'Qa' not in df.columns:
            raise ValueError("No 'Qa' column was found. Make sure you have the flow rates of the pipes under a column named 'Qa' ")
        
        if not K and equation in ['darcy', 'hazen']:
            if 'length' not in df.columns:
                raise ValueError("No 'length' column was found. Make sure you have the lengths of the pipes under a column named 'length' ")
            elif 'diameter' not in df.columns:
                raise ValueError("No 'diameter' column was found. Make sure you have the Diameter of the pipes under a column named 'diameter' ")
            elif equation == 'darcy' and 'f' not in df.columns:
                raise ValueError("No 'f' column was found. Make sure you have the friction factor(f) of the pipes under a column named 'f' ")
            elif equation == 'hazen' and 'C' not in df.columns:
                raise ValueError("No 'C' column was found. Make sure you have the Hazen Williams's coefficient(C) of the pipes under a column named 'C'")
        else:
            if 'K' not in df.columns:
                raise ValueError("No 'K' column was found. Make sure you have the K values of the pipes under a column named 'K' ")
        
        # Copy the assumed flow rate to a new column to be used for the calculation
        df['Qnew'] = df[f'Qa']
        df_list.append({'name': lp['name'], 'df': df})
        
    iterations = 1
    while True:
        q_list = []   
        
        # Function to calculate the head loss in each pipe
        def head_loss(a, b, c, d):
            # a=length, b=diameter, c=f/C/K, d=assumed rate
            
            hl = math.pow(abs(d), n)
            if not K and equation == 'darcy':
                const_A = 8*a*c; 
                const_B = 9.81*math.pow(3.143, 2)*math.pow(b, 5)
                hl = hl*(const_A/const_B)
            elif not K and equation == 'hazen':
                const_A = 10.67*a
                const_B = math.pow(c, n)*math.pow(b, 4.87)
                hl = hl*(const_A/const_B)
            else:
                hl = hl*c
            
            return hl if d > 0 else -hl
        
        for df in df_list:
            if not K and equation == 'darcy':
                df['df']['hl'] = df['df'].apply(lambda row: head_loss(row['length'], row['diameter'], row['f'], row['Qnew']), axis=1)
            elif not K and equation == 'hazen':
                df['df']['hl'] = df['df'].apply(lambda row: head_loss(row['length'], row['diameter'], row['C'], row['Qnew']), axis=1)
            else:
                df['df']['hl'] = df['df'].apply(lambda row: head_loss(1, 1, row['K'], row['Qnew']), axis=1)
            
            # Divide the head loss by the rate
            df['df']['hl/Qnew'] = abs(df['df']['hl']/(df['df']['Qnew']))
            
            # Find the sum
            hl_sum = sum(df['df']['hl'])
            hl_Qnew_sum = sum(df['df']['hl/Qnew'])
            
            # calculate the error, check if it's less than 0.2 and adjust the rate
            q = -hl_sum/(1.85*hl_Qnew_sum)
            q_list.append({'name': df['name'], 'q': q})
        
        common_loop_pipes = []
        if pipes_in_common_loops is not None:
            for com_pipe in pipes_in_common_loops:
                common_loop_pipes.append(com_pipe['name'])
                
        def correct_rate(row, df, df_q):
            Qnew = row['Qnew'] + df_q
            if pipes_in_common_loops and row['pipe'] in common_loop_pipes:
                common_p = None
                
                for com_p in pipes_in_common_loops:
                    if com_p['name'] == row['pipe']:
                        common_p = com_p
                
                if common_p:
                    common_pipe_loops = []
                    for loop_name in common_p['loops']:
                        if loop_name != df['name']:
                            common_pipe_loops.append(loop_name)
                            
                    if len(common_pipe_loops) != 0:
                        for loop_name in common_pipe_loops:
                            for q in q_list:
                                if loop_name == q['name']:
                                    Qnew = Qnew - q['q']
                
            return Qnew
        
        df_q_list = []
        for df_q in q_list:
            df_q_list.append(df_q['q'])
        
        if all(abs(item) < error for item in df_q_list):
            break
        
        for df, df_q in zip(df_list, q_list):     
            df['df']['Qnew'] = df['df'].apply(lambda row: correct_rate(row, df, df_q['q']), axis=1)
        
        if iterations >= max_iter:
            print("No convergence after 100 iterations. You can increase the number of iterations by modifying the 'max_iter' parameter in the function argument")
            break

        iterations +=1
    
    for df, df_q in zip(df_list, q_list):
        df['df']['Qnew'] = round(df['df']['Qnew'], 4)
        df['df']['hl'] = round(df['df']['hl'], 4)
        df['df']['hl/Qnew'] = round(df['df']['hl/Qnew'], 4)
        df_to_list = df['df'].values.tolist()
        df_to_list.insert(0, df['df'].columns.tolist())
        df['df'] = df_to_list
        df['props'] = df.pop('df')

    return df_list



