function [ distinct ] = kernel_num(mem,order )
% Numero di componenti distinti in un nucleo di Volterra
% di memoria mem e ordine order

distinct = factorial(mem)/factorial(mem-order)/factorial(order)+mem^order-factorial(mem)/factorial(mem-order);
end

