{
  description = "minimal independent set calculator and CNF minimizer";
  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixpkgs-unstable";
    arjun = {
      url = "github:itepastra/arjun/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    cryptominisat = {
      url = "github:itepastra/cryptominisat/add-flake";
      inputs.nixpkgs.follows = "nixpkgs";
    };
    sbva = {
      url = "github:itepastra/sbva/add-flake";
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
          zlib,
          cryptominisat,
          arjun,
          sbva,
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
            zlib
            cryptominisat
            arjun
            sbva
          ];
          postInstall = ''mv $out/include/approxmc/approxmc.h $out/include/approxmc.h'';
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
          };
        in
        {
          inherit approxmc;
          default = approxmc;
        }
      );
    };
}
