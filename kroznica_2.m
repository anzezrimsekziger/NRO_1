function kroznica_2()
    st = 10000;
    r = 1;
    
    [ocena, napaka] = area_pi(st, r);
    fprintf('Ocenjeno π: %f\nNapaka: %f\n', ocena, napaka);
    
    risi_kroznico_in_tocke(st, r);
end

function [ocena, napaka] = area_pi(st, r)
    xy = 2 * rand(2, st) - 1;
    razdalja = sum(xy.^2);
    znotrajKroga = sum(razdalja <= r^2);
    
    ocena = 4 * znotrajKroga / st;
    napaka = abs(ocena - pi);
end

function risi_kroznico_in_tocke(st, r)
    xy = 2 * rand(2, st) - 1;
    razdalja = sum(xy.^2);
    
    tockeN = razdalja <= r^2;
    tockeZ = razdalja > r^2;
    
    scatter(xy(1, tockeN), xy(2, tockeN), 50, 'b', 'filled');
    hold on;
    scatter(xy(1, tockeZ), xy(2, tockeZ), 50, 'r', 'x');
    
    theta = linspace(0, pi/2, 1000);
    x = r*cos(theta);
    y = r*sin(theta);
    plot(x, y);
    
    %axis equal;
    %title('Naključno generirane točke na krožnici z lokom');
    %xlabel('X-os');
    %ylabel('Y-os');
    %legend('Znotraj krožnice', 'Zunaj krožnice', 'Krožnica');
end

