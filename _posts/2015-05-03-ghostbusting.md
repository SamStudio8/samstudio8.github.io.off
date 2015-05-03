---
layout: post
title: "[meta] Ghostbusting"
---

Shortly after setting up this blog, I embedded Google Analytics tracking; primarily because I like numbers
but also in hope of discovering that at least one other person who isn't me or one my supervisors is interested
in my adventures. It's also great writing practice and ensures that I properly think through the things that
I am doing to avoid looking [wrong on the internet](https://xkcd.com/386/).

I was already in the habit of spamming links to my posts via various social networks so it wasn't a long
wait for the warm, fuzzy feeling of confirmation that people were actually reading my work. Or at the very
least accessing it and bouncing away.

However after a few days, I noticed several strange entries:

| Source                          | Sessions | % of Referrals | Bounce Rate | Pages / Session | Avg. Session Duration |
|---------------------------------|----------|----------------|-------------|-----------------|-----------------------|
| pornhub-forum.uni.me            | 90       | 37.50%         | 97.78%      | 1.02            | 00:00:05              |
| free-share-buttons.com          | 79       | 32.92%         | 8.86%       | 1.91            | 00:01:24              |
| site4.free-share-buttons.com    | 25       | 10.42%         | 0.00%       | 2.00            | 00:01:32              |
| site3.free-share-buttons.com    | 18       | 7.50%          | 0.00%       | 2.00            | 00:01:27              |
| site2.free-share-buttons.com    | 17       | 7.08%          | 0.00%       | 2.00            | 00:01:28              |
| forum.topic62206786.darodar.com | 7        | 2.92%          | 0.00%       | 3.00            | 00:00:00              |
| www.Get-Free-Traffic-Now.com    | 4        | 1.62%          | 100.00%     | 1.00            | 00:00:00              |

Curses. All my non-social referrals are **ghost referrals**! It seems that spam bots can cycle through guesses for Google Analytics
tracking IDs and remotely execute the tracking script and so hit without even needing to access the website.
