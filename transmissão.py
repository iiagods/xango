import numpy as np
import matplotlib.pyplot as plt

class ModulacaoEscalar:
    def __init__(self, e_star, mu_values, num_cycles=1):
        self.e_star = e_star
        self.mu_values = mu_values
        self.num_cycles = num_cycles
        self.t = None
        self.vs10_star = None
        self.vs20_star = None
        self.vs30_star = None
        self.i1 = None
        self.i2 = None
        self.i3 = None
        self.thd_values = []
        self.f_onda_p = 50
        self.f_onda_m = 5
        self.pi23 = 2 * np.pi / 3
        self.rq23 = np.sqrt(2 / 3)
        self.rq3 = np.sqrt(3)
        
    def t_entrada(self):
        self.vs1_star = np.sin(2 * np.pi * self.t)
        self.vs2_star = np.sin(2 * np.pi * self.t - 2 * np.pi / 3)
        self.vs3_star = np.sin(2 * np.pi * self.t + 2 * np.pi / 3)

    def calculate_vN0max_vN0mim(self):
        self.vN0max_star = self.e_star / 2 - np.maximum.reduce([self.vs1_star, self.vs2_star, self.vs3_star])
        self.vN0mim_star = -self.e_star / 2 - np.minimum.reduce([self.vs1_star, self.vs2_star, self.vs3_star])

    def calculates_vN0(self, mu):
        self.vN0_star = mu * self.vN0max_star + (1 - mu) * self.vN0mim_star

    def tm(self):
        self.vs10_star = self.vs1_star + self.vN0_star
        self.vs20_star = self.vs2_star + self.vN0_star
        self.vs30_star = self.vs3_star + self.vN0_star

    def c(self):
        self.i1 = np.sin(2 * np.pi * self.t)
        self.i2 = np.sin(2 * np.pi * self.t - 2 * np.pi / 3)
        self.i3 = np.sin(2 * np.pi * self.t + 2 * np.pi / 3)

    def calculate_thd(self): # Placeholder para cálculo real do THD
        thd = np.random.random() * 10  # Exemplo: THD fictício entre 0 e 10%
        return thd

    def simulate_mod(self, mu):
        self.t = np.linspace(0, self.num_cycles, 1000)  # Vetor de tempo
        self.t_entrada()  # Calcular tensões de entrada
        self.calculate_vN0max_vN0mim()  # Calcular vN0max* e vN0mim*
        self.calculates_vN0(mu)  # Calcular vN0*
        self.tm()  # Calcular tensões moduladas
        self.c()  # Calcular correntes

    def plot_mod(self, mu):
        plt.figure(figsize=(10, 6))
        
        plt.subplot(2, 1, 1)
        plt.plot(self.t, self.vs10_star, label=r'$v_{s10}^*$', color='black')
        plt.title(f'Tensão Modulada $v_{{s10}}^*$ (µ = {mu})')
        plt.ylabel('Tensão')
        plt.legend()
        
        plt.subplot(2, 1, 2)
        plt.plot(self.t, self.i1, label=r'$i_1$', color='blue')
        plt.plot(self.t, self.i2, label=r'$i_2$', color='red')
        plt.plot(self.t, self.i3, label=r'$i_3$', color='green')
        plt.title('Correntes')
        plt.ylabel('Corrente')
        plt.legend()
        
        plt.tight_layout()
        plt.show()

    def plot_thd_vs_mu(self):
        plt.figure(figsize=(8, 6))
        plt.plot(self.mu_values, self.thd_values, label='THD', color='black', marker='o')
        plt.title('THD em função de µ')
        plt.xlabel('µ')
        plt.ylabel('THD (%)')
        plt.legend()
        plt.grid(True)
        plt.show()

    def run(self):
        for mu in self.mu_values:
            self.simulate_mod(mu)
            thd = self.calculate_thd()
            self.thd_values.append(thd)
            self.plot_mod(mu)

def exemplo_modulacao_escalar():
    e_star = 1.0  # Tensão DC do inversor
    mu_values = [0.5, 1.0]  # Valores de µ
    simulation = ModulacaoEscalar(e_star, mu_values)
    simulation.run()

# Executa o exemplo
exemplo_modulacao_escalar()
