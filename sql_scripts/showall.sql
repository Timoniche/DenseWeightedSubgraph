SELECT DonorInfo.*, Densities.density
FROM DonorInfo INNER JOIN Densities
ON Donorinfo.info_id = Densities.info_id