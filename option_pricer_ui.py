from option_pricer_sub_tabs import * 

class OptionPricerApp:
    def __init__(self, root):
        self.root = root
        self.root.title('Option Pricer')

        # Create a style for tabs
        self.style = ttk.Style()
        self.style.theme_create('SleekStyle', parent='alt', settings={
            'TNotebook': {'configure': {'tabmargins': [2, 5, 2, 0]}},
            'TNotebook.Tab': {
            'configure': {'padding': [10, 5], 'background': '#2E8B57', 'foreground': 'white'},
            'map': {'background': [('selected', '#3CB371')], 
                    'expand': [('selected', [1, 1, 1, 0])]}}})
        self.style.theme_use('SleekStyle')

        # Create a notebook (tabs container)
        self.notebook = ttk.Notebook(root)
        self.notebook.pack(fill=tk.BOTH, expand=True)

        # Define functions for each tab
        self.tab_functions = [
            ('European Options', tab_black_scholes_options),
            ('Implied volatility Calculator', tab_iv_calculator),
            ('Geometric Asian Options', tab_geometric_asian_options_closed_form),
            ('Arithmetic Asian Options', tab_arithmetic_asian_options),
            ('Geometric Basket Options', tab_geometric_basket_options_closed_form),
            ('Arithmetic Basket Options', tab_arithmetic_basket_options),
            ('KIKO Put Options', tab_kiko_options),
            ('American Options', tab_american_options_binomial_tree),
        ]

        # Create tabs            
        for name, function in self.tab_functions:
            tab = ttk.Frame(self.notebook)
            self.notebook.add(tab, text=name)
            function(tab)

def main():
    root = tk.Tk()
    app = OptionPricerApp(root)
    root.mainloop()

if __name__ == '__main__':
    main()
