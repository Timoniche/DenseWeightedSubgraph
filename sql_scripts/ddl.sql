CREATE TABLE Proximity
(
    function_id INT NOT NULL PRIMARY KEY,
    function_code VARCHAR(1000) NOT NULL
);

CREATE TABLE DonorInfo
(
    info_id SERIAL NOT NULL,
    donor_id VARCHAR(100) NOT NULL,
    chr INT NOT NULL,
    function_id INT NOT NULL,
    UNIQUE (donor_id, chr, function_id),
    PRIMARY KEY (info_id),
    FOREIGN KEY (function_id) REFERENCES Proximity (function_id)
);

CREATE TABLE Clusters
(
    info_id INTEGER NOT NULL,
    bp INT NOT NULL,
    PRIMARY KEY (info_id),
    FOREIGN KEY (info_id) REFERENCES DonorInfo (info_id)
);

CREATE TABLE Periphery
(
    info_id INTEGER NOT NULL,
    bp INT NOT NULL,
    PRIMARY KEY (info_id),
    FOREIGN KEY (info_id) REFERENCES DonorInfo (info_id)
);

CREATE TABLE Densities
(
    info_id INTEGER NOT NULL,
    density REAL NOT NULL,
    PRIMARY KEY (info_id),
    FOREIGN KEY (info_id) REFERENCES DonorInfo (info_id)
);