.onLoad <- function(libname, pkgname)
{
	library.dynam(package="networkReliability", chname="networkReliability", lib.loc = .libPaths())
}
