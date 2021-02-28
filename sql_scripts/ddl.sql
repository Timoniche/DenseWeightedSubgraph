CREATE TABLE Proximity
(
    FunctionId INT NOT NULL PRIMARY KEY,
    FunctionCode VARCHAR(1000) NOT NULL
);

CREATE TABLE Clusters
(
    DonorId VARCHAR(100) NOT NULL,
    ChrN INT NOT NULL,
    FunctionId INT NOT NULL,
    Bp INT NOT NULL,
    PRIMARY KEY (DonorId, ChrN, FunctionId),
    FOREIGN KEY (FunctionId) REFERENCES Proximity (FunctionId)
);

CREATE TABLE Periphery
(
    DonorId VARCHAR(100) NOT NULL,
    ChrN INT NOT NULL,
    FunctionId INT NOT NULL,
    Bp INT NOT NULL,
    PRIMARY KEY (DonorId, ChrN, FunctionId),
    FOREIGN KEY (FunctionId) REFERENCES Proximity (FunctionId)
);

CREATE TABLE Densities
(
    DonorId VARCHAR(100) NOT NULL,
    ChrN INT NOT NULL,
    FunctionId INT NOT NULL,
    Density REAL NOT NULL,
    PRIMARY KEY (DonorId, ChrN, FunctionId),
    FOREIGN KEY (FunctionId) REFERENCES Proximity (FunctionId)
);