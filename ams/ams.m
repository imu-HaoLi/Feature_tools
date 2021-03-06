function cepDpDD = ams(sub_gf)

nChan = 64;

%sub_gf = reshape(sub_gf, 1, length(sub_gf));
sub_gf = sub_gf';

nFrame = floor(length(sub_gf)/160)-1;

ns_ams = extract_AMS_perChan(sub_gf,nChan,16e3);

AMS_points = ns_ams';
AMS_points = AMS_points(1:nFrame,:);

AMS_points = AMS_points.';

%AMS_points = [zeros(size(AMS_points,1),1) AMS_points];

cepDpDD = AMS_points;

end
