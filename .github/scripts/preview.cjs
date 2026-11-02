const fs = require('node:fs');
const path = require('node:path');

// This script is loaded from the default branch by the privileged publisher.
// Artifact metadata is data, never executable JavaScript or shell source.
module.exports = async ({github, context, core, artifactRoot = 'public'}) => {
  core.setOutput('publish', 'false');
  const run = context.payload.workflow_run;
  if (run.event !== 'pull_request' || run.conclusion !== 'success') return;
  const metadata = JSON.parse(fs.readFileSync(path.join(artifactRoot, 'pr.json'), 'utf8'));
  if (!Number.isSafeInteger(metadata.number) || metadata.number <= 0) {
    throw new Error('Invalid PR number in preview artifact');
  }
  const {owner, repo} = context.repo;
  const {data: pull} = await github.rest.pulls.get({owner, repo, pull_number: metadata.number});
  const repository = `${owner}/${repo}`;
  // Production-origin previews are only for code from trusted in-repository branches.
  // Fork builds remain downloadable artifacts, even if their verification was approved.
  if (pull.state !== 'open' || pull.base.repo.full_name !== repository ||
      pull.head.repo?.full_name !== repository || pull.head.sha !== run.head_sha ||
      run.head_repository?.full_name !== repository) {
    core.info('Skipping closed, stale, or fork preview');
    return;
  }
  const site = path.join(artifactRoot, 'site');
  if (!fs.lstatSync(site).isDirectory()) throw new Error('Preview site must be a directory');
  const inspect = directory => {
    for (const entry of fs.readdirSync(directory, {withFileTypes: true})) {
      if (entry.name.toLowerCase() === '.git' || (!entry.isFile() && !entry.isDirectory())) {
        throw new Error('Preview contains Git metadata, a link, or a special file');
      }
      if (entry.isDirectory()) inspect(path.join(directory, entry.name));
    }
  };
  inspect(site);
  core.setOutput('number', String(pull.number));
  core.setOutput('url', `https://lib.cp-algorithms.com/${pull.number}/`);
  core.setOutput('publish', 'true');
};
