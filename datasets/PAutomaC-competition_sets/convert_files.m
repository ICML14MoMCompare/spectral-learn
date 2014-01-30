problems = {'1','6','7','14','27','29','33','39','42','43','45','26'};

for i=1:12
    interpret_pautomac(strcat(problems(i),'.pautomac_model.txt'), strcat(problems(i),'.wfa.text'));
end