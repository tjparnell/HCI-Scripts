package RT::Interface::Email::HCIWhiteListFilter;

our $VERSION = '1.0';

=head1 NAME

RT::Interface::Email::HCIWhiteListFilter - only accept new tickets via email from whitelist

=head1 SYNOPSIS

Doesn't accept new emails to the support queue from outside reasonable domains.
In other words, simple, crude spam filter.

=head1 INSTALL

	# edit etc/RT_SiteConfig.pm

	Set(@Plugins,(qw/
		RT::Interface::Email::HCIWhiteListFilter
	/));
	Set(@MailPlugins, (qw/
		Auth::MailFrom
		HCIWhiteListFilter
	/));
	Set(%Plugin_HCIWhiteListFilter, (
		"pass"	=> [ qr/\@\hci\.utah\.edu$/i, qr/\.utah\.edu$/i ], 
		# "message" => "Error: Please submit requests from appropriate email domain.",
   ));

The pass value is an array of C<qr> regex expressions for acceptable domains.

If you want an error message emailed back to the sender, put in a message.

=head1 CREDIT

Basic design and implementation modeled after RT::Interface::Email::RequiredHeaders
by Alister West.

=head1 AUTHOR

 Timothy J. Parnell, PhD
 Dept of Oncological Sciences
 Huntsman Cancer Institute
 University of Utah
 Salt Lake City, UT, 84112

=head1 LICENCE AND COPYRIGHT

Copyright 2015, Timothy Parnell

This package is free software; you can redistribute it and/or modify
it under the terms of the Artistic License 2.0.  

=cut

use 5.008;
use warnings;
use strict;

use RT::Interface::Email qw(ParseSenderAddressFromHead);

=head1 GetCurrentUser - RT MailPlugin Callback

	Returns: ($CurrentUser, $auth_level) - not-triggered passthough inputs.
	Returns: ($CurrentUser, -1 )		 - halt further processing and send rejection notice.

	See RT::Interface::Email::GetAuthenticationLevel for more on $auth_level.

=cut

sub GetCurrentUser {
	my %args = (
		Message		  => undef,
		RawMessageRef => undef,
		CurrentUser	  => undef,
		AuthLevel	  => undef,
		Action		  => undef,
		Ticket		  => undef,
		Queue		  => undef,
		@_
	);

	# Default return values - if not triggering this plugin.
	my @ret = ( $args{CurrentUser}, $args{AuthLevel} );

	my %config = RT->Config->Get('Plugin_HCIWhiteListFilter');
	my $passDomains = $config{pass};
	my $message = $config{message};

	$RT::Logger->debug("HCIWhiteListFilter debugging ...");

	# If required was not supplied - skip this plugin.
	if (!$passDomains || !@$passDomains) {
		$RT::Logger->debug( " .. no 'pass' domains were set - SKIP");
		return @ret;
	}

	# we only ever filter 'new' tickets.
	if ($args{'Ticket'}->id) {
		$RT::Logger->debug( " .. ticket correspondence - SKIP");
		return @ret;
	}

	# check email sender domain
	my $head = $args{'Message'}->head;
	my ($address, $name, $error) = ParseSenderAddressFromHead($head);
	my $ok = 0;
	foreach my $domain (@$passDomains) {
		if ($address =~ /$domain/) {
			$ok = 1;
			last;
		}
	}
	
	# verify 
	if (not $ok) {
		
		# notify sender
		if ($message) {
			my $ErrorsTo = RT::Interface::Email::ParseErrorsToAddressFromHead($head);
			RT::Interface::Email::MailError(
				To			=> $ErrorsTo,
				Subject		=> "Permission denied : " . $head->get('Subject'),
				Explanation => $message,
# 					MIMEObj		=> $args{'Message'}
			);
			# don't contribute to email flood by attaching the original email
		}

		# halt further email processing to block creation of a ticket.
		$RT::Logger->info("HCIWhiteListFilter: [error] email from $address - HALT");
		return ( $args{CurrentUser}, -1 );
	}

	$RT::Logger->info("HCIWhiteListFilter: email domains ok");
	return @ret;
}


1;
