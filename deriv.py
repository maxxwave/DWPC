import numpy as np
import matplotlib.pyplot as plt

# Read data from file
file_path = 'x_t3.txt'
data = np.loadtxt(file_path)

# Calculate the derivative
derivative = np.gradient(data)

# Plot the original data and its derivative
plt.figure(figsize=(8, 6))

plt.subplot(2, 1, 1)
plt.plot(data, label='Original Data')
plt.title('Original Data')
plt.legend()

plt.subplot(2, 1, 2)
plt.plot(derivative, label='Derivative')
plt.title('Derivative of Data')
plt.legend()

plt.tight_layout()
plt.show()

lya=[]



def lyapunov(x_input,lya):
    sum_lya=0.0
    for i in range(0, len(x_input),1):
        lya.append(np.log(abs(x_input[i])))
        sum_lya += np.log(abs(x_input[i]))
        print("lyapunov coeff is..", lya[i])
    return sum_lya

print("Lyapunov_coeff is:....", lyapunov(derivative,lya)/len(derivative) )

plt.plot(lya, label='Lyapunov')
plt.title('lyap')
plt.show()


