---
layout: post
title: "[meta] Ghostbusting"
excerpt: "Shortly after setting up this blog, I embedded Google Analytics tracking; primarily because I like numbers. Sadly someone appears to be messing with my numbers."
---

Shortly after setting up this blog, I embedded Google Analytics tracking; primarily because I like numbers
but also in hope of discovering that at least one other person who isn't me or one my supervisors is interested
in my adventures. It's also great writing practice and gives me the chance to properly think through the things that
I am doing to avoid looking [wrong on the internet](https://xkcd.com/386/).

I was already in the habit of spamming links to my posts via various social networks so it wasn't a long
wait for the warm, fuzzy feeling of confirmation that people were actually reading my work. Or at the very
least clicking on it.

However after a few days, I noticed several strange entries amongst my lovely numbers[^0]:

| Source                          | Sessions | % of Referrals | Bounce Rate | Pages / Session | Avg. Session Duration |
|---------------------------------|----------|----------------|-------------|-----------------|-----------------------|
| pornhub-forum.uni\[dot\]me            | 90       | 37.50%         | 97.78%      | 1.02            | 00:00:05              |
| free-share-buttons\[dot\]com          | 79       | 32.92%         | 8.86%       | 1.91            | 00:01:24              |
| site4.free-share-buttons\[dot\]com    | 25       | 10.42%         | 0.00%       | 2.00            | 00:01:32              |
| site3.free-share-buttons\[dot\]com    | 18       | 7.50%          | 0.00%       | 2.00            | 00:01:27              |
| site2.free-share-buttons\[dot\]com    | 17       | 7.08%          | 0.00%       | 2.00            | 00:01:28              |
| forum.topic62206786.darodar\[dot\]com | 7        | 2.92%          | 0.00%       | 3.00            | 00:00:00              |
| Get-Free-Traffic-Now\[dot\]com        | 4        | 1.62%          | 100.00%     | 1.00            | 00:00:00              |

Curses. All my non-social referrals are **ghost referrals**! Disreputable publishers use spambots liberally to remotely
execute Google Analytics tracking scripts[^1] to appear to be providing a stream of referrals to your
website. Though, I'm unsure of the aim of this apparent **data pollution[^2]** attack. Beyond a poor attempt
at driving confused hostmaters to the sources to increase organic traffic I don't really see
what the benefit to the executor is[^3]? Perhaps the sites attempt to install malware on or capture more
valuable information from unsuspecting visitor's machines.

Oddly, page specific metrics (such as landing/exit pages) are polluted too. Bots copy the hostname of the
referral source to the page name of the false hit, giving hostmasters the impression something more worrying
is afoot. It's easy to forget falsification of data is a potential possibility, especially when one is not
responsible for collection and management of the data.

None of this is particularly important or bothersome, unless like me, you like numbers and numbers that are
wrong are upsetting. So how can normality be restored? These spambots target indiscriminately and remotely,
leaving them unaware of the actual target and thus with no option but to spoof the `hostname` of the hit
(or leave the field unset) which should in fact match that of the website under attack.

A [helpful blogpost](http://www.ohow.co/stop-adult-referral-spam-ga/#How_to_Block_these_Referrals_in_Google_Analytics) details how to set up a simple filter
in your Google Analytics control panel to remove future spurious data by ignoring hits which fail to provide
a valid expected hostname. Hooray.

* * *

# tl;dr
* Spambots performed a seemingly pointless **data pollution** attack on my Google Analytics records.
* One should always be as suspicious of data as possible, especially if it was collected by somebody else.

[^0]: I haven't censored the source column as it may be potentially useful to others having the same problem.

[^1]: Presumably by guessing or spidering Google Analytics tracking IDs.

[^2]: Electing to use "pollution" over "poison" here as the result is less directly toxic and more of
    confusing annoyance.
    
[^3]: Although by mentioning the domains here perhaps I've done exactly what they wanted...
