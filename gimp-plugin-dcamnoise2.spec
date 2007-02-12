
%define		plugin_name	dcamnoise2

Summary:	Plugin for removing noise introduced by digital cameras
Summary(pl.UTF-8):   Wtyczka usuwająca szumy wprowadzane przez aparaty cyfrowe
Name:		gimp-plugin-%{plugin_name}
Version:	0.63
Release:	1
License:	GPL
Group:		X11/Applications/Graphics
# http://registry.gimp.org/file/dcamnoise2-0.63.c?action=download&id=6329
Source0:	%{plugin_name}-%{version}.c
URL:		http://registry.gimp.org/plugin?id=5610
BuildRequires:	gimp-devel >= 1:2.0
Requires:	gimp >= 1:2.0
BuildRoot:	%{tmpdir}/%{name}-%{version}-root-%(id -u -n)

%define		_plugindir	%(gimptool --gimpplugindir)/plug-ins

%description
Plugin for removing noise introduced by digital cameras. Very
effective. It is present in image menu, under Filters/Enhance.

%description -l pl.UTF-8
Wtyczka usuwająca szumy wprowadzane przez aparaty cyfrowe. Bardzo
skuteczna. Jest dostępna w menu obrazka, w podmenu Filtry/Uwydatnianie
(ang. Filters/Enhance).

%prep
%setup -q -c -T

cp "%{SOURCE0}" "%{plugin_name}-%{version}.c"

%build
gimptool --build "%{plugin_name}-%{version}.c"

%install
rm -rf $RPM_BUILD_ROOT

install -d $RPM_BUILD_ROOT%{_plugindir}
install "%{plugin_name}-%{version}" $RPM_BUILD_ROOT%{_plugindir}

%clean
rm -rf $RPM_BUILD_ROOT

%files
%defattr(644,root,root,755)
%attr(755,root,root) %{_plugindir}/*
