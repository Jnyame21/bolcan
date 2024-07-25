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
    Str: A message indicating that the result did not converge.
    List: A list of dictionaries, each containing the loop's name and a list of data for each pipe, along with the assumed and corrected rates.
    
        
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
    results = True
    if equation not in ['darcy', 'hazen']:
        raise ValueError("The equation must either be 'hazen' for Hazen Williams or 'darcy' for Darcy Weisbach")
    
    n = 2 if equation == 'darcy' else 1.85
    
    df_list = []
    for lp in loops:
        header = lp['props'][0]
        data = lp['props'][1:]
        df = pd.DataFrame(data, columns=header)
        df['Qa'] = df['Qa'].astype(float)
        
        if 'pipe' not in df.columns:
            raise ValueError("No 'pipe' column was found. Make sure you have the names of the pipes under a column named 'pipe' ")
        elif 'Qa' not in df.columns:
            raise ValueError("No 'Qa' column was found. Make sure you have the flow rates of the pipes under a column named 'Qa' ")
        
        if not K:
            if 'length' not in df.columns:
                raise ValueError("No 'length' column was found. Make sure you have the lengths of the pipes under a column named 'length' ")
            elif 'diameter' not in df.columns:
                raise ValueError("No 'diameter' column was found. Make sure you have the Diameter of the pipes under a column named 'diameter' ")
            elif equation == 'darcy' and 'f' not in df.columns:
                raise ValueError("No 'f' column was found. Make sure you have the friction factor(f) of the pipes under a column named 'f' ")
            elif equation == 'hazen' and 'C' not in df.columns:
                raise ValueError("No 'C' column was found. Make sure you have the Hazen Williams's coefficient(C) of the pipes under a column named 'C'")
            else:
                df['length'] = df['length'].astype(float)
                df['diameter'] = df['diameter'].astype(float)
        else:
            if 'K' not in df.columns:
                raise ValueError("No 'K' column was found. Make sure you have the K values of the pipes under a column named 'K' ")
            else:
                df['K'] = df['K'].astype(float)
        
        if not K and equation == 'darcy':
            df['f'] = df['f'].astype(float)
        elif not K and equation == 'hazen':
            df['C'] = df['C'].astype(float)
            
        # Copy the assumed flow rate to a new column to be used for the calculation
        df['Qnew'] = df['Qa']
        df_list.append({'name': lp['name'], 'df': df})
        
    iterations = 1
    
    while True:
        q_list = []   
        
        # Function to calculate the head loss in each pipe
        def head_loss(a, b, c, d):
            # a=length, b=diameter, c=f/C/K, d=assumed rate
            a=float(a); b=float(b); c=float(c); d=float(d)
            
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
            
            # calculate the correction factor
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
        
        # Check if all the correction factors are less than the error tolerance
        if all(abs(item) < error for item in df_q_list):
            break
        
        for df, df_q in zip(df_list, q_list):     
            df['df']['Qnew'] = df['df'].apply(lambda row: correct_rate(row, df, df_q['q']), axis=1)
        
        if iterations >= max_iter:
            results = False
            break

        iterations +=1
    
    if results:
        for df, df_q in zip(df_list, q_list):
            df['df']['Qnew'] = round(df['df']['Qnew'], 4)
            df['df']['hl'] = round(df['df']['hl'], 4)
            df['df']['hl/Qnew'] = round(df['df']['hl/Qnew'], 4)
            df_to_list = df['df'].values.tolist()
            df_to_list.insert(0, df['df'].columns.tolist())
            df['df'] = df_to_list
            df['props'] = df.pop('df')

        return df_list
    else:
        return f"No convergence after {max_iter} iterations. You can increase the number of iterations by modifying the 'max_iter' parameter in the function argument"


def gas_velcity(Qb, Pb, Tb, D, T, z, P, unit:str):
    """
    Calculate the velocity of gas flowing through a pipeline
    
    Args:
        Qb(float): The base flow rate of the gas in m3/day or ft3/day(scf/day)
        D(float): The diameter of the pipe in mm or inches
        Pb(float): The base pressure of the gas in kPa or psia
        Tb(float): The base temperature of the gas in Kelvin(K) or Rankine(R)
        T(float): The temperature of the flowing gas in Kelvin(K) or Rankine(R)
        z(float): The gas compressibility factor, dimensionaless
        P(float): The pressure of the gas, kPa or psia
        
    Returns:
        float: The velocity of the gas in the pipeline
    
    Example
    ---------
    Qb = 250e6    # base flow rate, ft3/day (scf/day)
    Pb = 14.7     # base pressure, psia
    Tb = 520      # base temperature, Rankine
    D = 19        # diameter, inches
    T = 520       # temparature, Rankine
    z = 1         # z factor
    P = 1014.7      # pressure, psia
    unit = 'uscs' # unit
    
    result = gas_velcity(Qb=Qb, Pb=Pb, Tb=Tb, D=D, T=T, z=z, P=P, unit='uscs')
    print(result)
    
    Result
    21.289
    
    """
    if unit not in ['si', 'uscs']:
        raise ValueError("The unit parameter must either be 'si' for SI units or 'uscs' for US customary unit")
    
    v = (Qb/math.pow(D, 2))*(Pb/Tb)*((T*z)/P)
    if unit == 'uscs':
        const_A = 0.002122
        v = v*const_A
    else:
        const_A = 14.7349
        v = v*const_A
    
    return round(v, 3)


def zfactor_cnga(Pavg, G, Tf):
    """
    Calculate the compressiblity factor(z-factor) of a gas in a pipeline using the california natural gas association method(cnga)
    
    Args:
        P(float): Average pressure, psig
        G(float): Specific gravity, dimensionaless
        Tf(float): The gas flowing temperature, Rankine(R)
        
    Returns:
        float: The compressibility factor(z-factor) of the gas
    
    Example
    ---------
    G = 0.6
    Tf = 520
    Pavg = 1000
    result = zfactor_cnga(Pavg=Pavg, G=G, Tf=Tf)
    print(result)
    
    Result
    0.858
    
    """
    const_A = 1 + ((Pavg*344400*math.pow(10, 1.785*G))/math.pow(Tf, 3.825))
    return round(1/const_A, 3)


def reynold_number(Q, Pb, Tb, D, G, u, unit):
    """
    Calculate the Reynold number of gas flowing through a pipeline
    
    Args:
        Q(float): Flow rate, m3/day or ft3/day(scf/day)
        Pb(float): Base pressure, kPa or psia
        Tb(float): Base temperature, Kelvin(K) or Rankine(R)
        D(float): Pipe diameter, mm or inches
        G(float): Gas specify gravity, dimensionless
        u(float): Gas viscosity, lb/ft-s or poise
        unit(str): Unit, si or uscs
        
    Returns:
        float: The Reynold number
    
    Example
    ---------
    Q = 250e6     # flow rate, ft3/day (scf/day)
    Pb = 14.7     # base pressure, psia
    Tb = 520      # base temperature, Rankine
    D = 19        # diameter, inches
    G = 0.6       # gas gravity
    u = 0.00008   # viscosity
    unit = 'uscs' # unit
    
    result = reynold_number(Q=Q, Pb=Pb, Tb=Tb, D=D, G=G, u=u, unit='uscs')
    print(result)
    
    Result
    1332931.4
    
    """
    if unit not in ['si', 'uscs']:
        raise ValueError("The unit parameter must either be 'si' for SI units or 'uscs' for US customary unit")
    
    v = (Pb/Tb)*((Q*G)/(u*D))
    if unit == 'uscs':
        const_A = 0.0004778
        v = v*const_A
    else:
        const_A = 0.5134
        v = v*const_A
    
    return round(v, 1)



