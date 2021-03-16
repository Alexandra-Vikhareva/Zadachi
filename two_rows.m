function B = two_rows(A)
    if (size(A,1)==2) &  (ndims(A)==2)
        B = A
    else c = size(A)
         B = zeros(c)
         disp('Массив должен быть двумерным и иметь две строки!')
    end
end