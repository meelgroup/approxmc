{
  description = "Approximate model counter";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    arjun = {
      url = "github:meelgroup/arjun/aig-shared-ptr-better-fix";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:msoos/cryptominisat/working-on-synth";
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
  };
  outputs =
    {
      self,
      nixpkgs,
      arjun,
      cryptominisat,
      sbva,
      evalmaxsat,
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
          fetchFromGitHub,
          gmp,
          mpfr,
          zlib,
          cryptominisat,
          arjun,
          sbva,
          evalmaxsat,
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
            ];
          };
          nativeBuildInputs = [
            cmake
            autoPatchelfHook
          ];
          buildInputs = [
            gmp
            mpfr
            zlib
            cryptominisat
            arjun
            sbva
            evalmaxsat
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
          };
        in
        {
          inherit approxmc;
          default = approxmc;
        }
      );
    };
}
