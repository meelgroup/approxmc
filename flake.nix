{
  description = "Approximate model counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    arjun = {
      url = "github:meelgroup/arjun/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:msoos/cryptominisat/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    sbva = {
      url = "github:meelgroup/sbva/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    evalmaxsat = {
      url = "github:meelgroup/EvalMaxSAT/master";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    treedecomp = {
      url = "github:meelgroup/treedecomp/main";
      inputs.nixpkgs.follows = "nixpkgs";
    };
  };
  outputs =
    {
      self,
      nixpkgs,
      arjun,
      cryptominisat,
      sbva,
      evalmaxsat,
      treedecomp,
    }:
    let
      inherit (nixpkgs) lib;
      systems = lib.intersectLists lib.systems.flakeExposed lib.platforms.linux;
      forAllSystems = lib.genAttrs systems;
      nixpkgsFor = forAllSystems (system: nixpkgs.legacyPackages.${system});
      fs = lib.fileset;
      approxmc-package =
        {
          stdenv,
          cmake,
          autoPatchelfHook,
          pkg-config,
          fetchFromGitHub,
          gmp,
          mpfr,
          zlib,
          cryptominisat,
          arjun,
          sbva,
          evalmaxsat,
          treedecomp,
        }:
        stdenv.mkDerivation {
          name = "approxmc";
          src = fs.toSource {
            root = ./.;
            fileset = fs.unions [
              ./CMakeLists.txt
              ./cmake
              ./src
              ./approxmcConfig.cmake.in
              ./pyproject.toml
            ];
          };
          nativeBuildInputs = [
            cmake
            autoPatchelfHook
            pkg-config
          ];
          cmakeFlags = [
            "-Dcryptominisat5_DIR=${cryptominisat}/lib/cmake/cryptominisat5"
            "-Darjun_DIR=${arjun}/lib/cmake/arjun"
          ];
          buildInputs = [
            gmp
            mpfr
            zlib
            cryptominisat
            arjun
            sbva
            evalmaxsat
            treedecomp
          ];
        };
    in
    {
      packages = forAllSystems (
        system:
        let
          approxmc = nixpkgsFor.${system}.callPackage approxmc-package {
            arjun = arjun.packages.${system}.arjun;
            cryptominisat = cryptominisat.packages.${system}.cryptominisat;
            sbva = sbva.packages.${system}.sbva;
            evalmaxsat = evalmaxsat.packages.${system}.evalmaxsat;
            treedecomp = treedecomp.packages.${system}.treedecomp;
          };
        in
        {
          inherit approxmc;
          default = approxmc;
        }
      );
    };
}
