function x = posroot(Roots)

% This subfunction extracts the positive real roots from
% those obtained in the call to MATLAB’s 'roots' function.
% If there is more than one positive root, the first one is chosen

%Roots: the vector of roots of a polynomial

%Construct the vector of positive real roots:
posroots = Roots(find(Roots>0 & ~imag(Roots)));
npositive = length(posroots);

%Exit if no positive roots exist:
if npositive == 0
  return
end
%...If there is more than one positive root, output the
% roots to the command window and prompt the user to
% select which one to use:
if npositive == 1
  x = posroots;
else
  x = posroots(1);
end

end