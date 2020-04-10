%% Function Test

% Set up a bunch of epsilon values
low = mp('10e7');
high = low*100;
inc = 990000;
smallSeq = low:inc:high;
eps = arrayfun(@(x) 1/x, smallSeq);

%% First Jacobi block... k = 3
k = 3;
J3 = mp([0 1 0; 0 0 1; 0 0 0]);

% Get corresponding radii
radEpsJ3 = mp(zeros(1, length(eps)));
for ind = 1:length(eps)
    radEpsJ3(ind) = findradiusMP(J3,eps(ind));
end

knownTerms = [1 1];
[p1j3, c1j3] = findNextPower(eps, radEpsJ3, k, knownTerms);
% Returns p1 = 1.666
%         c1 = 1

knownTerms = [1 1; 1 5/3];
[p2j3, c2j3] = findNextPower(eps, radEpsJ3, k, knownTerms);
% Returns p2 = 2.333
%         c2 = 1.666

knownTerms = [1 1; 1 5/3; 5/3 7/3];
[p3j3, c3j3] = findNextPower(eps, radEpsJ3, k, knownTerms);
% Returns p3 = 2.96 ~ 3 (which would be the next one in the series)
%         c3 = 0.178



%% Moving on... k=4
k = 4;
J4 = mp([0 1 0 0; 0 0 1 0; 0 0 0 1; 0 0 0 0]);

radEpsJ4 = mp(zeros(1, length(eps)));
for ind = 1:length(eps)
    radEpsJ4(ind) = findradiusMP(J4,eps(ind));
end

knownTerms = [1 1];
[p1j4, c1j4] = findNextPower(eps, radEpsJ4, k, knownTerms);
% p1 --> 1.5
% c1 --> 1

knownTerms = [1 1; 1 1.5];
[p2j4, c2j4] = findNextPower(eps, radEpsJ4, k, knownTerms);
% p2 --> 2
% c2 --> 1.5

knownTerms = [1 1; 1 1.5; 1.5 2];
[p3j4, c3j4] = findNextPower(eps, radEpsJ4, k, knownTerms);
% p3 --> 2.5
% c3 --> 2.625

knownTerms = [1 1; 1 1.5; 1.5 2; 2.625 2.5];
[p4j4, c4j4] = findNextPower(eps, radEpsJ4, k, knownTerms);


%% Moving on... k=5
k = 5;
J5 = mp([0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1; 0 0 0 0 0]);

radEpsJ5 = mp(zeros(1, length(eps)));
for ind = 1:length(eps)
    radEpsJ5(ind) = findradiusMP(J5,eps(ind));
end

knownTerms = [1 1];
[p1j5, c1j5] = findNextPower(eps, radEpsJ5, k, knownTerms);
% p1 --> 1.4
% c1 --> 1

knownTerms = [1 1; 1 1.4];
[p2j5, c2j5] = findNextPower(eps, radEpsJ5, k, knownTerms);
% p2 --> 1.8
% c2 --> 1.4

knownTerms = [1 1; 1 1.4; 1.4 1.8];
[p3j5, c3j5] = findNextPower(eps, radEpsJ5, k, knownTerms);
% p3 --> 2.2
% c3 --> 2.24

k=5;
knownTerms = [1 1; 1 1.4; 1.4 1.8; 2.24 2.2];
[p4j5, c4j5] = findNextPower(eps, radEpsJ5, k, knownTerms);


%% Moving on... k=6
k = 6;
J6 = mp([0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1; 0 0 0 0 0 0]);

radEpsJ6 = mp(zeros(1, length(eps)));
for ind = 1:length(eps)
    radEpsJ6(ind) = findradiusMP(J6,eps(ind));
end

knownTerms = [1 1];
[p1j6, c1j6] = findNextPower(eps, radEpsJ6, k, knownTerms);
% p1 --> 1.333
% c1 --> 1

knownTerms = [1 1; 1 8/6];
[p2j6, c2j6] = findNextPower(eps, radEpsJ6, k, knownTerms);


knownTerms = [1 1; 1 8/6; 8/6 10/6];
[p3j6, c3j6] = findNextPower(eps, radEpsJ6, k, knownTerms);

k=6;
knownTerms = [1 1; 1 8/6; 8/6 10/6; 2 2];
[p4j6, c4j6] = findNextPower(eps, radEpsJ6, k, knownTerms);
