function [vector] = string_to_double_vector(string)

element1 = '';
element2 = '';
element3 = '';

for i=2:length(string)-1
    if (string(i)~=' ')
        element1 = [element1 +string(i)];
    else
        j = i;
        break;
    end
end
for i=j+1:length(string)-1
    if (string(i)~=' ')
        element2 = [element2 +string(i)];
    else
        k = i;
        break;
    end
end
for i=k+1:length(string)-1
    if (string(i)~=' ')
        element3 = [element3 +string(i)];
    else
        break;
    end
end
vector(1) = str2double(element1);
vector(2) = str2double(element2);
vector(3) = str2double(element3);
end

