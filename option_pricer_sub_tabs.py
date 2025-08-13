import tkinter as tk
from tkinter import ttk
from option_pricer_functions import *

# helper function to create entry boxes and labels
def create_entry_grid(frame, labels):
    entries = {}
    for i, (label_text, var_name) in enumerate(labels):
        tk.Label(frame, text=label_text).grid(row=i, column=0, sticky='w')
        entry = tk.Entry(frame)
        entry.grid(row=i, column=1)
        entries[var_name] = entry
    return entries
    
def tab_black_scholes_options(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma = float(sigma_entry.get())
            r = float(r_entry.get())
            q = float(q_entry.get())
            OptionType = option_type_var.get()[0].upper()  # Get the first letter of the selected option type
            
            # Calculate option price
            option_price = black_scholes_options(S, K, T, sigma, r, q, OptionType)
            
            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f'Option Price: {option_price:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')
    
    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility (sigma):', 'sigma_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Repo Rate (q):', 'q_entry'),
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_entry = entries['sigma_entry']
    r_entry = entries['r_entry']
    q_entry = entries['q_entry']

    tk.Label(input_frame, text='Option Type:').grid(row=6, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'Call', 'Put')
    option_type_menu.grid(row=6, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=2, width=40, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)
    
    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_iv_calculator(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            r = float(r_entry.get())
            q = float(q_entry.get())
            option_premium = float(option_premium_entry.get())
            OptionType = option_type_var.get()[0].upper()

            # Calculate iv
            iv = iv_calculator(option_premium, S, K, T, r, q, OptionType)

            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f'Implied volatility: {iv:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Repo Rate (q):', 'q_entry'),
        ('Option Premium:', 'option_premium_entry')
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    r_entry = entries['r_entry']
    q_entry = entries['q_entry']
    option_premium_entry = entries['option_premium_entry']
    
    tk.Label(input_frame, text='Option Type:').grid(row=6, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=6, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=2, width=40, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_american_options_binomial_tree(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma = float(sigma_entry.get())
            r = float(r_entry.get())
            OptionType = option_type_var.get()[0].upper()
            n = int(N_entry.get())

            # Calculate option price
            option_price = american_options_binomial_tree(S, K, T, sigma, r, OptionType, n)
            
            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f'American Option Price: {option_price:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)
    
    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility (sigma):', 'sigma_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Number of Steps (N):', 'N_entry')  
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_entry = entries['sigma_entry']
    r_entry = entries['r_entry']
    N_entry = entries['N_entry']

    tk.Label(input_frame, text='Option Type:').grid(row=6, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=6, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=2, width=40, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_geometric_asian_options_closed_form(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma = float(sigma_entry.get())
            r = float(r_entry.get())
            OptionType = option_type_var.get()[0].upper()
            n = int(n_entry.get())

            # Calculate option price
            option_price = geometric_asian_options_closed_form(S, K, T, sigma, r, OptionType, n)

            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f'Geometric Asian Option Price: {option_price:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility (sigma):', 'sigma_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Number of Obervation times (n):', 'n_entry')  
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_entry = entries['sigma_entry']
    r_entry = entries['r_entry']
    n_entry = entries['n_entry']

    tk.Label(input_frame, text='Option Type:').grid(row=6, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=6, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=2, width=40, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_arithmetic_asian_options(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma = float(sigma_entry.get())
            r = float(r_entry.get())
            n = int(n_entry.get())
            OptionType = option_type_var.get()[0].upper()
            m = int(m_entry.get())
            use_control_variate = True if control_variate_var.get() == 'yes' else False

            # Calculate option price and 95% confidence interval
            res = arithmetic_asian_options(S, K, T, sigma, r, OptionType, n, monto_carlo_num_paths=m, use_control_variate=use_control_variate)

            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END,
                               f'Arithmetic Asian Option Price: {res[0]:.4f} \n'
                               f'95% confidence interval: [{res[1][0]:.4f}, {res[1][1]:.4f}]')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)
    
    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility (sigma):', 'sigma_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Number of Obervation times (n):', 'n_entry'),
        ('Number of Paths to Simulate (m):', 'm_entry')
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_entry = entries['sigma_entry']
    r_entry = entries['r_entry']
    n_entry = entries['n_entry']
    m_entry = entries['m_entry']
    
    tk.Label(input_frame, text='Option Type:').grid(row=7, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=7, column=1)

    tk.Label(input_frame, text='Use Control Variate:').grid(row=8, column=0, sticky='w')
    control_variate_var = tk.StringVar()
    control_variate_var.set('yes')
    control_variate_menu = tk.OptionMenu(input_frame, control_variate_var, 'yes', 'no')
    control_variate_menu.grid(row=8, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=4, width=50, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_geometric_basket_options_closed_form(tab):
    def calculate():
        try:
            S_1 = float(S1_entry.get())
            S_2 = float(S2_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma_1 = float(sigma_1_entry.get())
            sigma_2 = float(sigma_2_entry.get())
            correlation = float(correlation_entry.get())
            r = float(r_entry.get())
            OptionType = option_type_var.get()[0].upper()

            # Calculate option price
            option_price = geometric_basket_options_closed_form(S_1, S_2, K, T, sigma_1, sigma_2, correlation, r, OptionType)

            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, f'Geometric Basket Option Price: {option_price:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price of Asset 1 (S1):', 'S1_entry'),
        ('Spot Price of Asset 2 (S2):', 'S2_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility of Asset 1 (sigma1):', 'sigma_1_entry'),
        ('Volatility of Asset 2 (sigma2):', 'sigma_2_entry'),
        ('Correlation (rho):', 'correlation_entry'),
        ('Risk-free Rate (r):', 'r_entry')
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S1_entry = entries['S1_entry']
    S2_entry = entries['S2_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_1_entry = entries['sigma_1_entry']
    sigma_2_entry = entries['sigma_2_entry']
    correlation_entry = entries['correlation_entry']
    r_entry = entries['r_entry']

    tk.Label(input_frame, text='Option Type:').grid(row=8, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=8, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=2, width=40, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_arithmetic_basket_options(tab):
    def calculate():
        try:
            S_1 = float(S1_entry.get())
            S_2 = float(S2_entry.get())
            K = float(K_entry.get())
            T = float(T_entry.get())
            sigma_1 = float(sigma_1_entry.get())
            sigma_2 = float(sigma_2_entry.get())
            correlation = float(correlation_entry.get())
            r = float(r_entry.get())
            OptionType = option_type_var.get()[0].upper()
            m = int(m_entry.get())
            use_control_variate = True if control_variate_var.get() == 'yes' else False

            # Calculate option price
            res = arithmetic_basket_options(S_1, S_2, K, T, sigma_1, sigma_2, correlation, r, OptionType, monto_carlo_num_paths=m, use_control_variate=use_control_variate)
                
            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END,
                               f'Arithmetic Basket Option Price: {res[0]:.4f} \n'
                               f'95% confidence interval: [{res[1][0]:.4f}, {res[1][1]:.4f}]')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price of Asset 1 (S1):', 'S1_entry'),
        ('Spot Price of Asset 2 (S2):', 'S2_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility of Asset 1 (sigma1):', 'sigma_1_entry'),
        ('Volatility of Asset 2 (sigma2):', 'sigma_2_entry'),
        ('Correlation (rho):', 'correlation_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Number of Paths to Simulate (m):', 'm_entry')
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S1_entry = entries['S1_entry']
    S2_entry = entries['S2_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_1_entry = entries['sigma_1_entry']
    sigma_2_entry = entries['sigma_2_entry']
    correlation_entry = entries['correlation_entry']
    r_entry = entries['r_entry']
    m_entry = entries['m_entry']
    
    tk.Label(input_frame, text='Option Type:').grid(row=9, column=0, sticky='w')
    option_type_var = tk.StringVar()
    option_type_var.set('call')  # Default option type
    option_type_menu = tk.OptionMenu(input_frame, option_type_var, 'call', 'put')
    option_type_menu.grid(row=9, column=1)

    tk.Label(input_frame, text='Use Control Variate:').grid(row=10, column=0, sticky='w')
    control_variate_var = tk.StringVar()
    control_variate_var.set('yes')
    control_variate_menu = tk.OptionMenu(input_frame, control_variate_var, 'yes', 'no')
    control_variate_menu.grid(row=10, column=1)

    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=4, width=50, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)

def tab_kiko_options(tab):
    def calculate():
        try:
            S = float(S_entry.get())
            sigma = float(sigma_entry.get())
            r = float(r_entry.get())
            T = float(T_entry.get())
            K = float(K_entry.get())
            barrier_lower = float(barrier_lower_entry.get())
            barrier_upper = float(barrier_upper_entry.get())
            n = int(n_entry.get())
            rebate = float(rebate_entry.get())

            # Calculate option price
            res = kiko_options(S, K, T, sigma, r, barrier_lower, barrier_upper, rebate, n)
    
            # Update output text
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END,
                               f'KIKO Put Option Price: {res[0]:.4f} \n'
                               f'Delta: {res[1]:.4f}')
        except ValueError:
            output_text.delete(1.0, tk.END)
            output_text.insert(tk.END, 'Please check your inputs.')

    # Create input labels and entry boxes
    input_frame = ttk.Frame(tab)
    input_frame.pack(padx=20, pady=20)

    entries = create_entry_grid(input_frame, [
        ('Spot Price (S):', 'S_entry'),
        ('Strike (K):', 'K_entry'),
        ('Time to Maturity (T):', 'T_entry'),
        ('Volatility (sigma):', 'sigma_entry'),
        ('Risk-free Rate (r):', 'r_entry'),
        ('Lower Barrier (L):', 'barrier_lower_entry'),
        ('Upper Barrier (U):', 'barrier_upper_entry'),
        ('Cash Rebate (R):', 'rebate_entry'),
        ('Number of Obervation times (n):', 'n_entry')  
    ])

    # Access entries using their variable names, e.g., S_entry = entries['S_entry']
    S_entry = entries['S_entry']
    K_entry = entries['K_entry']
    T_entry = entries['T_entry']
    sigma_entry = entries['sigma_entry']
    r_entry = entries['r_entry']
    barrier_lower_entry = entries['barrier_lower_entry']
    barrier_upper_entry = entries['barrier_upper_entry']
    rebate_entry = entries['rebate_entry']
    n_entry = entries['n_entry']
    
    # Button to calculate option price
    calculate_button = ttk.Button(tab, text='Calculate', command=calculate, style='Border.TButton')
    calculate_button.pack(pady=10)

    # Output text box
    output_text = tk.Text(tab, height=4, width=50, wrap=tk.WORD)
    output_text.pack(padx=20, pady=10)

    # Add a line border to the Calculate Option button with modern styling
    tab_style = ttk.Style()
    tab_style.configure('Border.TButton', borderwidth=1, relief='solid', background='#f0f0f0', font=('Helvetica', 12), padding=5)